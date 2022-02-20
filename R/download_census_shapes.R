#' Download and Prepare Map Shapes
#'
#' Download cartographic boundary files from the [U.S. Census Bureau](https://www.census.gov),
#' and reformat them to simplified GeoJSON files.
#'
#' @param dir Path to a directory in which to save the reformatted files. If these are to be
#' used in a published site, they should be in the \code{"docs"} directory of the site directory. If not specified,
#' the shapes are downloaded and returned, but not saved.
#' @param fips State Federal Information Processing Standard (FIPS) code, United States Postal Service (USPS)
#' code, or name; see [state.txt](https://www2.census.gov/geo/docs/reference/state.txt).
#' @param entity Entity name (e.g., \code{"county"}, \code{"tract"}, or \code{"bg"}); see
#' [2019_file_name_def.pdf](https://www2.census.gov/geo/tiger/GENZ2019/2019_file_name_def.pdf)
#' (or appropriate year's file). States and counties are only available at the national level, so if these are
#' requested for a particular state, the national files will be downloaded and subset.
#' @param name Name for the GeoJSON file (without extension) to be written.
#' @param year Year of the shapes, 2010 and after.
#' @param resolution Resolution of the shapes; one of \code{"500k"} (default), \code{"5m"}, or \code{"20m"}.
#' @param strip_features Logical; if \code{TRUE}, will strip all features other than IDs.
#' @param simplify A function to simplify shapes with, such as \code{rmapshaper::ms_simplify}. This
#' function should take the \code{sf} \code{data.frame} as its first argument, and return that object
#' with its \code{geometry} replaced.
#' @param ... Passes additional arguments to \code{simplify}; e.g., for \code{rmapshaper::ms_simplify},
#' you might want to specify \code{keep} and set \code{keep_shapes} to \code{TRUE}.
#' @param force Logical; if \code{TRUE}, will force a re-download, and overwrite any existing files.
#' @examples
#' \dontrun{
#' # download Virginia counties shapes
#' shapes <- download_census_shapes(
#'   "docs/shapes", "va", "county",
#'   name = "counties",
#'   simplify = rmapshaper::ms_simplify, keep_shapes = TRUE
#' )
#' }
#' @return An \code{sf} \code{data.frame} with a geometry column containing shapes, and other features in columns.
#' @export

download_census_shapes <- function(dir = NULL, fips = "us", entity = "state", name = NULL, year = 2019,
                                   resolution = "500k", strip_features = FALSE, simplify = NULL, ..., force = FALSE) {
  us_fips <- list(
    name = c(
      "united states", "alabama", "alaska", "arizona", "arkansas", "california", "colorado", "connecticut",
      "delaware", "district of columbia", "florida", "georgia", "hawaii", "idaho", "illinois", "indiana", "iowa",
      "kansas", "kentucky", "louisiana", "maine", "maryland", "massachusetts", "michigan", "minnesota", "mississippi",
      "missouri", "montana", "nebraska", "nevada", "new hampshire", "new jersey", "new mexico", "new york",
      "north carolina", "north dakota", "ohio", "oklahoma", "oregon", "pennsylvania", "rhode island", "south carolina",
      "south dakota", "tennessee", "texas", "utah", "vermont", "virginia", "washington", "west virginia", "wisconsin",
      "wyoming", "american samoa", "guam", "northern mariana islands", "puerto rico", "minor outlying islands",
      "virgin islands"
    ),
    post = c(
      "us", "al", "ak", "az", "ar", "ca", "co", "ct", "de", "dc", "fl", "ga", "hi", "id", "il", "in", "ia", "ks", "ky",
      "la", "me", "md", "ma", "mi", "mn", "ms", "mo", "mt", "ne", "nv", "nh", "nj", "nm", "ny", "nc", "nd", "oh", "ok",
      "or", "pa", "ri", "sc", "sd", "tn", "tx", "ut", "vt", "va", "wa", "wv", "wi", "wy", "as", "gu", "mp", "pr", "um",
      "vi"
    ),
    fips = c(
      "us", 1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
      33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56, 60, 66, 69, 72, 74, 78
    )
  )
  entities <- c(
    "aiannh", "anrc", "bg", "cbsa", "cd116", "cnecta", "concity", "county", "county_within_cd116",
    "county_within_ua", "cousub", "csa", "division", "elsd", "nation", "necta", "place", "puma10", "region", "scsd",
    "sldl", "sldu", "state", "subbarrio", "tract", "ua10", "unsd", "zcta510"
  )
  fips <- gsub("^0|[^a-z0-9\\s]", "", tolower(fips[[1]]))
  fid <- which(us_fips$fips == fips)
  if (!length(fid)) {
    fid <- if (nchar(fips) == 2) {
      which(us_fips$post == fips)
    } else {
      pmatch(fips, us_fips$name)
    }
    if (!length(fid) || is.na(fid)) cli_abort("failed to recognize {.arg fips}")
  }
  fips <- us_fips$fips[fid]
  if (nchar(fips) == 1) fips <- paste0("0", fips)
  if (nchar(year) < 4) year <- as.numeric(paste0("20", year))
  if (year < 2010) cli_abort("only years 2010 and after are available from this function")
  entity <- tolower(entity)
  entity <- match.arg(entity, entities)
  subset <- NULL
  if (fips != "us" && (entity == "county" || entity == "state")) {
    subset <- fips
    fips <- "us"
  }
  resolution <- tolower(resolution)
  resolution <- match.arg(resolution, c("500k", "5m", "20m"))
  file <- paste0("cb_", year, "_", fips, "_", entity, "_", resolution)
  url <- paste0("https://www2.census.gov/geo/tiger/GENZ", year, "/shp/", file, ".zip")
  temp <- tempdir()
  temp_file <- paste0(temp, "/", file, ".zip")
  out_file <- paste0(dir, "/", if (!is.null(name)) name else file, ".geojson")
  if (force && file.exists(out_file)) unlink(out_file)
  if (!file.exists(out_file)) {
    if (!force && file.exists(paste0(temp, "/", file))) temp_file <- paste0(temp, "/", file)
    if (force || !file.exists(temp_file)) {
      download.file(url, temp_file)
      unzip(temp_file, exdir = paste0(temp, "/", file))
      file.remove(temp_file)
      temp_file <- paste0(temp, "/", file)
    }
    shapes <- st_read(temp_file, quiet = TRUE)
    if (!is.null(subset)) shapes <- shapes[shapes$STATEFP == subset, ]
    if (strip_features) shapes <- shapes[, "GEOID", drop = FALSE]
    if (is.function(simplify)) shapes <- simplify(shapes, ...)
    if (!is.null(dir)) {
      if (!file.exists(dir)) dir.create(dir, FALSE, TRUE)
      st_write(shapes, out_file)
    }
  } else {
    shapes <- st_read(out_file, quiet = TRUE)
  }
  invisible(shapes)
}
