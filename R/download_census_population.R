#' Download U.S. Census Population Data
#'
#' Download and load U.S. census American Community Survey 5-year summary files at the block group level,
#' optionally including home-work commute statistics.
#'
#' @param dir Directory in which to save the file(s).
#' @param state Name, abbreviation, or FIPS code of the state.
#' @param year 4 digit year, between 2009 and the most recent year.
#' @param include_margins Logical; if \code{TRUE}, will save and include estimate margins or error, which is a
#' \code{data.frame} with the same dimensions as the estimates \code{data.frame}.
#' @param include_commutes Logical; if \code{TRUE}, will download the Longitudinal Employer-Household Dynamics (LEHD)
#' Origin-Destination Employment Statistics (LODES) data for the state and block groups within selected tracts.
#' @param counties,tracts,blockgroups A vector of counties, tracts, or block group GEOIDs within the specified state to
#' filter block groups for (block groups within the specified counties or tracts, or matching the specified block group
#' GEOID). Only one can be specified (lower levels overwrite higher levels). Defaults to all block groups.
#' @param overwrite Logical; if \code{TRUE}, will remove any existing files before downloading and saving new versions.
#' @param verbose Logical; if \code{FALSE}, will not print status messages.
#' @examples
#' \dontrun{
#' download_census_population(".", "va")
#' }
#' @return A list with at least an \code{estimates} entry, and optionally \code{margins} and/or \code{commutes}
#' entries. The \code{extimates} and \code{margins} entries are \code{data.frame}s with block groups in rows,
#' population variables in columns, and person counts in cells. The \code{commutes} entry is a sparse matrix with
#' home block groups in rows, work block groups in columns, and job counts in cells (all jobs, summed across blocks).
#' @export

download_census_population <- function(dir, state, year = 2019, include_margins = FALSE, include_commutes = FALSE,
                                       counties = NULL, tracts = NULL, blockgroups = NULL, overwrite = FALSE, verbose = TRUE) {
  us_fips <- list(
    name = c(
      "UnitedStates", "Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut",
      "Delaware", "DistrictOfColumbia", "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa",
      "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi",
      "Missouri", "Montana", "Nebraska", "Nevada", "NewHampshire", "NewJersey", "NewMexico", "NewYork",
      "NorthCarolina", "NorthDakota", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "RhodeIsland", "SouthCarolina",
      "CouthDakota", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", "WestWirginia", "Wisconsin",
      "Wyoming", "PuertoRico"
    ),
    post = c(
      "us", "al", "ak", "az", "ar", "ca", "co", "ct", "de", "dc", "fl", "ga", "hi", "id", "il", "in", "ia", "ks", "ky",
      "la", "me", "md", "ma", "mi", "mn", "ms", "mo", "mt", "ne", "nv", "nh", "nj", "nm", "ny", "nc", "nd", "oh", "ok",
      "or", "pa", "ri", "sc", "sd", "tn", "tx", "ut", "vt", "va", "wa", "wv", "wi", "wy", "pr"
    ),
    fips = c(
      "us", 1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
      33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56, 72
    )
  )
  state <- gsub("^0|[^a-z0-9]", "", tolower(state[[1]]))
  fid <- which(us_fips$fips == state)
  if (!length(fid)) {
    fid <- if (nchar(state) == 2) {
      which(us_fips$post == state)
    } else {
      pmatch(state, tolower(us_fips$name))
    }
    if (!length(fid) || is.na(fid)) cli_abort("failed to recognize {.arg state}")
  }
  state <- us_fips$name[fid]
  state_post <- us_fips$post[fid]
  display_state <- gsub("([a-z])([A-Z])", "\\1 \\2", state)

  url <- paste0(
    "https://www2.census.gov/programs-surveys/acs/summary_file/", year,
    "/data/5_year_seq_by_state/", state, "/Tracts_Block_Groups_Only/"
  )
  files_zip <- paste0(year, "5", state_post, "000", 1:3, "000.zip")
  if (!is.null(dir)) {
    dir <- tdir <- paste0(normalizePath(dir, "/", FALSE), "/")
    if (!dir.exists(dir)) dir.create(dir, FALSE, TRUE)
  } else {
    tdir <- paste0(tempdir(), "/")
  }
  folder <- paste0(tdir, state, "_Tracts_Block_Groups_Only/")
  files <- paste0(folder, "e", year, "5", state_post, "000", 1:3, "000.txt")
  final_files <- paste0(if (is.null(dir)) tdir else dir, state_post, "_", c(
    "population",
    if (include_margins) "population_margins"
  ), ".csv")
  if (overwrite) unlink(c(paste0(tdir, files_zip), files, final_files))
  if (!all(file.exists(final_files)) && !all(file.exists(paste0(tdir, files_zip)) | file.exists(files))) {
    if (verbose) cli_alert_info(paste("downloading", display_state, "American Community Survey summary files"))
    for (i in seq_along(files)) {
      if (!file.exists(files[i])) {
        download.file(paste0(url, files_zip[i]), paste0(tdir, files_zip[i]))
      }
    }
  }
  if (any(file.exists(paste0(tdir, files_zip)))) {
    if (verbose) cli_alert_info(paste("decompressing", display_state, "American Community Survey summary files"))
    dir.create(folder, FALSE)
    for (f in files_zip[file.exists(paste0(tdir, files_zip))]) {
      system2("unzip", c("-qd", shQuote(folder), shQuote(paste0(tdir, f))))
      unlink(paste0(tdir, f))
    }
  }

  if (!all(file.exists(final_files))) {
    if (verbose) cli_alert_info(paste("reformatting", display_state, "American Community Survey summary files"))

    # column names
    headers <- read.csv(paste0(
      "https://www2.census.gov/programs-surveys/acs/summary_file/", year,
      "/documentation/user_tools/ACS_5yr_Seq_Table_Number_Lookup.txt"
    ), nrows = 661)
    title_id <- which(headers$Total.Cells.in.Table != "")
    names <- sub("^Universe:\\s*", "", headers$Table.Title)
    sub_id <- which(grepl(":", names, fixed = TRUE))

    headers$name <- sub("\\.$", "", sub("_\\w+_Total", "_Total", gsub("._", "_", paste0(
      rep(gsub("[;:(),\\s]+", ".", names[title_id], perl = TRUE), c(title_id[-1], nrow(headers)) - title_id),
      "_",
      rep(gsub("[;:(),\\s]+", ".", names[sub_id], perl = TRUE), c(sub_id[-1], nrow(headers)) - sub_id),
      "_",
      gsub("[;:(),\\s]+", ".", sub(" years", "", names, fixed = TRUE), perl = TRUE)
    ), fixed = TRUE)))
    start_id <- which(!is.na(headers$Start.Position))
    offset <- headers$Line.Number - 1 + rep(
      headers[start_id, "Start.Position"],
      c(start_id[-1], nrow(headers) + 1) - start_id
    )

    ind_start <- which(!is.na(offset == min(offset, na.rm = TRUE)) & offset == min(offset, na.rm = TRUE))
    ind_end <- c(ind_start[-1], length(offset)) - ind_start
    inds <- lapply(seq_along(ind_start), function(i) {
      is <- seq(ind_start[i], ind_end[i])
      is[!is.na(offset[is]) & offset[is] %% 1 == 0]
    })

    subs <- as.numeric(sub("\\D.*$", "", headers$Total.Cells.in.Table))
    ind_start <- which(
      !is.na(headers$Start.Position) &
        headers$Start.Position == min(headers$Start.Position, na.rm = TRUE)
    )
    ind_end <- c(ind_start[-1] - 1, length(subs))
    inds <- lapply(seq_along(ind_start), function(i) {
      is <- seq(ind_start[i], ind_end[i])
      is[!is.na(headers[is, "Line.Number"]) & headers[is, "Line.Number"] %% 1 == 0]
    })

    # row geoids
    rows <- read.csv(paste0(
      "https://www2.census.gov/programs-surveys/acs/summary_file/", year,
      "/data/5_year_seq_by_state/", state, "/Tracts_Block_Groups_Only/g", year, "5", state_post, ".csv"
    ), header = FALSE)
    rows$GEOID <- sub("\\D.*$", "", paste0(
      formatC(rows$V10, width = 2, flag = 0),
      formatC(rows$V11, width = 3, flag = 0),
      formatC(rows$V14, width = 6, flag = 0),
      rows$V15
    ))

    # estimates
    estimates <- do.call(cbind, lapply(seq_along(files), function(i) {
      tab <- read.csv(files[i], header = FALSE)
      colnames(tab)[c(6, offset[inds[[i]]])] <- c("GEOID", headers[inds[[i]], "name"])
      tab$GEOID <- rows[tab$GEOID, "GEOID"]
      tab[nchar(tab$GEOID) == 12, -seq(1, if (i == 1) 5 else 6)]
    }))
    estimates[estimates == "."] <- NA

    # margins
    if (include_margins) {
      files_m <- paste0(folder, "m", year, "5", state_post, "000", 1:3, "000.txt")
      margins <- do.call(cbind, lapply(seq_along(files_m), function(i) {
        tab <- read.csv(files_m[i], header = FALSE)
        colnames(tab)[c(6, offset[inds[[i]]])] <- c("GEOID", headers[inds[[i]], "name"])
        tab$GEOID <- rows[tab$GEOID, "GEOID"]
        tab[nchar(tab$GEOID) == 12, -seq(1, if (i == 1) 5 else 6)]
      }))
      margins[margins == "."] <- NA
    }

    su <- vapply(estimates, function(v) all(is.na(v)), TRUE)
    if (any(su)) {
      estimates <- estimates[, !su, drop = FALSE]
      colnames(estimates) <- sub("\\.1$", "", colnames(estimates))
      if (include_margins) {
        margins <- margins[, !su, drop = FALSE]
        colnames(margins) <- sub("\\.1$", "", colnames(margins))
      }
    }

    su <- !is.na(estimates[, 2]) & estimates[, 2] == 0
    if (any(su)) {
      estimates <- estimates[!su, ]
      if (include_margins) margins <- margins[!su, ]
    }

    subsetter <- !c(is.null(blockgroups), is.null(tracts), is.null(counties))
    if (any(subsetter)) {
      su <- if (subsetter[1]) {
        estimates$GEOID %in% blockgroups
      } else if (subsetter[2]) {
        substr(estimates$GEOID, 1, 11) %in% tracts
      } else if (subsetter[3]) {
        substr(estimates$GEOID, 1, 5) %in% counties
      }
      if (any(su) && !all(su)) {
        estimates <- estimates[su, ]
        if (include_margins) margins <- margins[su, ]
      }
    }
  } else {
    if (verbose) cli_alert_info(paste("loading existing", display_state, "population data"))
    estimates <- read.csv(final_files[1])
    if (include_margins) margins <- read.csv(final_files[2])
  }


  # commutes
  if (include_commutes) {
    file <- paste0(tdir, state_post, "_commutes.csv")
    filegz <- paste0(file, ".gz")
    if (overwrite) unlink(c(file, filegz))
    if (!any(file.exists(c(file, filegz)))) {
      if (verbose) cli_alert_info(paste("downloading", display_state, "Origin-Destination Employment Statistics file"))
      download.file(paste0(
        "https://lehd.ces.census.gov/data/lodes/LODES7/", state_post,
        "/od/", state_post, "_od_main_JT00_", year, ".csv.gz"
      ), filegz)
    }
    if (fresh <- file.exists(filegz)) system2("gzip", c("-df", filegz))

    check <- file.exists(file)
    if (check && fresh && Sys.which("openssl")[[1]] != "") {
      hashes <- paste(readLines(paste0(
        "https://lehd.ces.census.gov/data/lodes/LODES7/", state_post, "/lodes_", state_post, ".sha256sum"
      )), collapse = "")
      hash_loc <- regexec(paste0("_od_main_JT00_", year), hashes, fixed = TRUE)[[1]]
      check <- substr(hashes, hash_loc - 68, hash_loc - 5) == strsplit(
        system2("openssl", c("dgst", "-sha256", shQuote(file)), TRUE), " ",
        fixed = TRUE
      )[[1]][2]
      if (!check) cli_warn("integrity check failed for {.file {file}}")
    }
    if (check) {
      if (verbose) cli_alert_info(paste0("loading ", if (!fresh) "existing ", display_state, " commutes data"))
      commutes <- read.csv(file)
      if (!is.null(commutes$w_geocode)) {
        fresh <- TRUE
        if (verbose) {
          cli_alert_info(paste(
            "reformatting", display_state, "Origin-Destination Employment Statistics file"
          ))
        }
        commutes$w_geocode <- substr(commutes$w_geocode, 1, 12)
        commutes$h_geocode <- substr(commutes$h_geocode, 1, 12)
        commutes <- commutes[commutes$w_geocode %in% estimates$GEOID & commutes$h_geocode %in% estimates$GEOID, ]
        ids <- structure(seq_along(estimates$GEOID), names = estimates$GEOID)
        edges <- tapply(commutes$S000, paste0(commutes$w_geocode, commutes$h_geocode), sum)
        commutes <- sparseMatrix(
          i = ids[substr(names(edges), 13, 24)],
          j = ids[substr(names(edges), 1, 12)],
          x = edges,
          dimnames = list(estimates$GEOID, estimates$GEOID)
        )
      } else if (colnames(commutes)[1] == "GEOID") {
        rownames(commutes) <- colnames(commutes)[-1] <- commutes[, 1]
        commutes <- as(as.matrix(commutes[, -1]), "dgCMatrix")
      }
    }
  }

  if (!is.null(dir)) {
    files <- paste0(dir, state_post, "_", c(
      "population",
      if (include_margins) "population_margins",
      if (include_commutes && check) "commutes"
    ), ".csv")
    if (!all(file.exists(files)) || (include_commutes && fresh && check)) {
      if (verbose) cli_alert_info(paste("writing", display_state, "files"))
      ck <- c(!file.exists(files[1]), !file.exists(files[2]) && include_margins, include_commutes && fresh && check)
      if (ck[1]) write.csv(estimates, files[1], row.names = FALSE)
      if (ck[2]) write.csv(margins, files[2], row.names = FALSE)
      if (ck[3]) write.csv(cbind(GEOID = rownames(commutes), as.matrix(commutes)), files[3], row.names = FALSE)
      if (verbose) {
        cli_bullets(c(
          v = paste0("wrote file", if (sum(ck) == 1) ":" else "s:"),
          structure(paste0("{.file ", files[ck], "}"), names = rep("*", sum(ck)))
        ))
      }
    }
  }

  invisible(list(estimates = estimates, margins = margins, commutes = commutes))
}
