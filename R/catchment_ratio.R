#' Calculate Floating Catchment Area Ratios
#'
#' Calculate provider (supply) to consumer (demand) ratios within floating catchment areas.
#'
#' @param consumers Number of consumers (demand); either a vector with consumer amounts (such as population counts),
#' or a matrix-like object with a column of IDs (located by \code{consumers_id}) and a column of amounts (located by
#' \code{consumers_value}).
#' @param providers Number of providers (supply); either a vector with provider amounts (such as number of doctors),
#' or a matrix-like object with a column of IDs (located by \code{providers_id}) and a column of amounts (located by
#' \code{providers_value}).
#' @param cost A matrix-like object of cost associated with each pair of \code{consumers} and \code{providers}
#' (such as distance or travel times), with \code{consumers} in rows and \code{providers} in columns.
#' \code{cost}'s dimensions should be \code{c(length(consumers), length(providers))} if \code{consumers} and
#' \code{providers} are vectors, or will be aligned by name if available (as vector names or in \code{_id} columns).
#' If \code{NULL}, coordinate information will be looked for in \code{consumers} and \code{providers} (based on
#' \code{consumers_location} and \code{providers_location}), from which to calculate Euclidean distances.
#' @param weight Means of defining catchment areas and their topology / friction. The simplest is a single number
#' representing a maximum distance between \code{consumers} and \code{providers} (2-step floating catchment area).
#' An enhancement of this is a list of vectors with two values each: the first is a distance, and the second a weight
#' to associate with that distance (e.g., \code{list(c(10, 1), c(20, .5))}, which will give consumers within a
#' \code{cost} of 10 full weight, and those within a \code{cost} of 20 half weight; enhanced 2-step floating catchment
#' area). If a character, refers to a weighting function (kernel density 2-step floating catchment area; in order
#' from most gradual to steepest between costs of \code{1} and \code{6}):
#' \tabular{ll}{
#'   \code{gaussian} (\code{ga}) \tab \code{exp(-cost^2 / (2 * scale ^ 2))}\cr
#'   \code{d*} (name of a density function; e.g., \code{"dnorm"}) \tab
#'     \code{weight(cost, 0, scale)}
#'   \cr
#'   \code{p*} (name of a distribution function; e.g., \code{"pnorm"}) \tab \code{weight(cost, 0, scale)}\cr
#'   \code{gravity} / \code{normal} (\code{gr} or \code{n}) \tab \code{sqrt(1 / cost^scale)}\cr
#'   \code{logarithmic} (\code{loga}) \tab \code{1 / (1 + log(cost, scale))}\cr
#'   \code{logistic} (\code{l}) \tab \code{1 / (1 + exp(scale * cost))}\cr
#'   \code{exponential} (\code{e}) \tab \code{exp(-cost * scale)}\cr
#' }
#' If a function, this will be passed \code{cost} as its first argument -- its output should be a matrix
#' convertible to a sparse matrix, of the same dimensions as \code{cost}. If a matrix-like object,
#' this will be converted to a sparse matrix.
#' @param normalize_weight Logical; if \code{TRUE}, weight is row-normalized such that \code{consumers} weights
#' are spread across \code{providers} in range. This can help correct for the increased weight of \code{consumers}
#' when they are in range of multiple \code{providers}. Selection weights like this make the difference between 2-
#' and 3-step floating catchment areas.
#' @param scale Numeric scaling factor if \code{weight} is the name of a decay function.
#' @param max_cost Numeric limit on \code{cost}. This is the same as setting \code{weight} to a single value,
#' or specifying a list of steps as \code{weight} (where the most distant step is effectively \code{max_cost}),
#' although a single-value weight is exclusive (\code{cost < weight}) where steps are inclusive. This is most useful
#' when \code{weight} is a weighing function, where \code{max_cost} will trim the tail of the weight distribution.
#' @param return_type Determines the values that are returned: \code{"original"} (default) for \code{providers}
#' per \code{consumers} (e.g., how many, likely fractional, doctors are accessible by each person within each region),
#' \code{"region"} for number of \code{providers} per \code{consumers} entry (\code{consumers * original}; e.g.,
#' how many doctors are accessible within each region), or \code{"normalized"} for \code{original} divided by
#' \code{sum(region) / sum(consumers)}.
#' @param consumers_id,consumers_value,consumers_location,providers_id,providers_value,providers_location Column
#' names in \code{consumers} and/or \code{providers} to extract IDs, values, and location data (referring to a single
#' \code{sf} geometry column, or multiple columns with coordinates). These can also be used to directly enter
#' ID, value, and/or location vectors (or matrices for location coordinates).
#' @param print_type Logical; if \code{TRUE}, will print the type of floating catchment area that was calculated.
#' @examples
#' pop <- c(5, 10, 50)
#' doc <- c(50, 100)
#' travel_time <- matrix(c(5, 50, 25, 70, 40, 30), ncol = 2)
#'
#' # 2-step floating catchment area
#' catchment_ratio(pop, doc, travel_time, 30)
#'
#' # kernel density (Gaussian) 2-step floating catchment area
#' catchment_ratio(pop, doc, travel_time, "gaussian")
#'
#' # enhanced 2-step floating catchment area
#' catchment_ratio(pop, doc, travel_time, list(c(60, .22), c(40, .68), c(20, 1)))
#'
#' # 3-step floating catchment area
#' catchment_ratio(pop, doc, travel_time, list(c(60, .22), c(40, .68), c(20, 1)), TRUE)
#'
#' # visualized weight functions
#' if (require("splot", quietly = TRUE)) {
#'   cost <- 1:10
#'   scale <- 2
#'   splot(list(
#'     gaussian = exp(-cost^2 / (2 * scale^2)),
#'     dnorm = dnorm(cost, 0, scale),
#'     pnorm = pnorm(-cost, 0, scale),
#'     gravity = sqrt(1 / cost^scale),
#'     logarithmic = 1 / (1 + log(cost, scale)),
#'     logistic = 1 / (1 + exp(scale * cost)),
#'     exponential = exp(-cost * scale)
#'   ) ~ cost, title = "Decay Functions", laby = "Weight", labx = "Cost", lines = "con", note = FALSE)
#' }
#' @return A vector with an access score (determined by \code{return_type}) for each entry in \code{consumers}.
#' @references
#' Dai, D. (2010). Black residential segregation, disparities in spatial access to health care facilities, and
#' late-stage breast cancer diagnosis in metropolitan Detroit. \emph{Health & place, 16}, 1038-1052.
#' doi: \href{https://doi.org/10.1016/j.healthplace.2010.06.012}{10.1016/j.healthplace.2010.06.012}
#'
#' Lou, W. & Qi, Y. (2009). An enhanced two-step floating catchment area (E2SFCA) method for measuring spatial
#' accessibility to primary care physicians. \emph{Health & Place, 15}, 1100-1107.
#' doi: \href{https://doi.org/10.1016/j.healthplace.2009.06.002}{10.1016/j.healthplace.2009.06.002}
#'
#' Luo, W. & Wang, F. (2003). Measures of spatial accessibility to health care in a GIS environment: synthesis
#' and a case study in the Chicago region. \emph{Environment and Planning B: Planning and Design, 30}, 865-884.
#' doi: \href{https://doi.org/10.1068/b29120}{10.1068/b29120}
#'
#' Wan, N., Zou, B., & Sternberg, T. (2012). A three-step floating catchment area method for analyzing spatial
#' access to health services. \emph{International Journal of Geographical Information Science, 26}, 1073-1089.
#' doi: \href{https://doi.org/10.1080/13658816.2011.624987}{10.1080/13658816.2011.624987}
#' @export

catchment_ratio <- function(consumers = NULL, providers = NULL, cost = NULL, weight = NULL, normalize_weight = FALSE,
                            scale = 2, max_cost = NULL, return_type = "original", consumers_id = "geoid",
                            consumers_value = "count", consumers_location = "geometry", providers_id = "geoid",
                            providers_value = "count", providers_location = "geometry", print_type = FALSE) {
  type <- ""
  input_data <- c(is.null(dim(providers)), is.null(dim(consumers)))
  # getting provider and consumer value vectors
  pv <- if (input_data[1]) {
    if(is.null(providers)){
      if(is.numeric(providers_value)) providers_value else
        cli_abort("{.arg providers} is not specified, and {.arg providers_value} is non-numeric")
    }else as.numeric(providers)
  } else if (length(providers_value) == 1 && providers_value %in% colnames(providers)) {
    as.numeric(providers[, providers_value, drop = TRUE])
  } else {
    as.numeric(providers[, which(vapply(seq_len(ncol(providers)), function(i) {
      is.numeric(providers[, i, drop = TRUE])
    }, TRUE))[1], drop = TRUE])
  }
  if (!length(pv)) cli_abort("failed to recognize values in {.arg providers}")
  cv <- if (input_data[2]) {
    if(is.null(consumers)){
      if(is.numeric(consumers_value)) consumers_value else
        cli_abort("{.arg consumers} is not specified, and {.arg consumers_value} is non-numeric")
    }else as.numeric(consumers)
  } else if (length(consumers_value) == 1 && consumers_value %in% colnames(consumers)) {
    as.numeric(consumers[, consumers_value, drop = TRUE])
  } else {
    as.numeric(consumers[, which(vapply(seq_len(ncol(consumers)), function(i) {
      is.numeric(consumers[, i, drop = TRUE])
    }, TRUE))[1], drop = TRUE])
  }
  if (!length(cv)) cli_abort("failed to recognize values in {.arg consumers}")
  # getting provider and consumer ids
  pid <- if (input_data[1]) {
    if (!is.null(names(providers))) names(providers) else if(length(providers) == length(providers_id))
      providers_id else seq_along(pv)
  } else if (length(providers_id) == 1 && providers_id %in% colnames(providers))
    providers[, providers_id, drop = TRUE] else seq_along(pv)
  if (!length(pid)) cli_abort("failed to recognize IDs in {.arg providers}")
  cid <- if (input_data[2]) {
    if (!is.null(names(consumers))) names(consumers) else if(length(consumers) == length(consumers_id))
      consumers_id else seq_along(cv)
  } else if (length(consumers_id) == 1 && consumers_id %in% colnames(consumers))
    consumers[, consumers_id, drop = TRUE] else seq_along(cv)
  if (!length(cid)) cli_abort("failed to recognize IDs in {.arg consumers}")
  if (is.null(cost)) {
    if ((!is.character(providers_location) || all(providers_location %in% colnames(providers))) &&
        (!is.character(consumers_location) || all(consumers_location %in% colnames(consumers)))) {
      pcords <- if(is.character(providers_location)) providers[, providers_location, drop = TRUE] else
        providers_location
      if (any(grepl("^sf", class(pcords)))) {
        if (!any(grepl("POINT$", class(pcords)))) pcords <- st_centroid(pcords)
        pcords <- st_coordinates(pcords)
      }
      ccords <- if(is.character(consumers_location)) consumers[, consumers_location, drop = TRUE] else
        consumers_location
      if (any(grepl("^sf", class(ccords)))) {
        if (!any(grepl("POINT$", class(ccords)))) ccords <- st_centroid(ccords)
        ccords <- st_coordinates(ccords)
      }
      if (ncol(pcords) == ncol(ccords)) {
        cost <- 1 / lma_simets(ccords, pcords, metric = "euclidean", pairwise = TRUE) - 1
      }else cli_abort(
        "{.arg cost} is NULL, and failed to calculate it from provided locations (differing number of columns)"
      )
    }
    if (is.null(cost)) {
      cost <- sparseMatrix(
        {},
        {},
        dims = c(length(cv), length(pv)),
        dimnames = list(cid, pid)
      )
      cost[] <- TRUE
    }
  } else if (is.null(dim(cost))) cli_abort("{.arg cost} must be a matrix-like object")
  if (is.null(weight) && !is.null(max_cost)) {
    weight <- max_cost
    max_cost <- NULL
  }
  if (is.null(weight)) {
    w <- as(cost > 0, "lgCMatrix")
  } else if (is.null(dim(weight))) {
    if (is.numeric(weight)) {
      # single buffer value means a uniformly weighted catchment area (original)
      w <- as(cost > 0 & cost < weight[[1]], "lgCMatrix")
    } else if (is.list(weight)) {
      # list of steps for roughly graded weightings (enhanced)
      type <- "enchanced "
      weight <- weight[order(-vapply(weight, "[[", 1, 1))]
      w <- as((cost <= weight[[1]][1]) * weight[[1]][2], "dgCMatrix")
      for (s in weight[-1]) w[cost <= s[1]] <- s[2]
      w[cost <= 0] <- 0
    } else if (is.character(weight)) {
      # name of a weight function (kernel density)
      weight <- tolower(weight[[1]])
      type <- paste0("kernel density (", weight, ") ")
      if ("units" %in% class(cost)) {
        cost <- matrix(
          as.numeric(cost), nrow(cost), ncol(cost),
          dimnames = dimnames(cost)
        )
      }
      if (grepl("^(?:gr|n)", weight)) {
        # gravity / normal kernel
        w <- as(sqrt(1 / cost^scale), "dgCMatrix")
      } else if (grepl("^e", weight)) {
        # exponential kernel
        w <- as(exp(-cost * scale), "dgCMatrix")
      } else if (grepl("^loga", weight)) {
        # logarithmic kernel
        w <- as(1 / (1 + log(cost, scale)), "dgCMatrix")
      } else if (grepl("^l", weight)) {
        # logistic kernel
        w <- as(1 / (1 + exp(scale * cost)), "dgCMatrix")
      } else if (grepl("^ga", weight)) {
        # Gaussian kernel
        w <- as(exp(-cost^2 / (2 * scale^2)), "dgCMatrix")
      } else if (grepl("^d", weight) && exists(weight)) {
        # assumed to be a density function like dnorm
        w <- as(cost, "dgCMatrix")
        w@x <- do.call(weight, list(w@x, 0, scale))
      } else if (grepl("^p", weight) && exists(weight)) {
        # assumed to be a distribution function like pnorm
        w <- as(cost, "dgCMatrix")
        w@x <- do.call(weight, list(-w@x, 0, scale))
      } else {
        cli_abort("{.arg weight} not recognized")
      }
      w[cost <= 0 | !is.finite(w)] <- 0
    } else if (is.function(weight)) {
      type <- paste0("kernel density (custom) ")
      w <- as(weight(cost), "dgCMatrix")
      w[cost <= 0 | !is.finite(w)] <- 0
    }
  } else {
    w <- as(weight, "dgCMatrix")
  }
  if (!is.null(max_cost)) w[cost > max_cost] <- 0
  wr <- rowSums(w)
  if(!is.null(rownames(w)) && !all(cid == rownames(w)) && all(cid %in% rownames(w))) w <- w[cid,]
  if(!is.null(colnames(w)) && !all(pid == colnames(w)) && all(pid %in% colnames(w))) w <- w[, pid]
  if(nrow(w) != length(cv) && is.numeric(cid) && min(cid) > 0 && max(cid) <= nrow(w)) w <- w[cid,]
  if(ncol(w) != length(pv) && is.numeric(pid) && min(pid) > 0 && max(pid) <= ncol(w)) w <- w[, pid]
  if (!all(dim(w) == c(length(cv), length(pv))))
    cli_abort("failed to align weight matrix with consumer and provider values")
  if (normalize_weight) {
    # adjust weights by selection probability (3-step)
    type <- paste0(type, "3-step floating catchment area")
    if(!is.null(rownames(w)) && all(rownames(w) %in% names(wr))) wr <- wr[rownames(w)]
    wr[wr == 0] <- 1
    w <- w * sweep(w, 1, wr, "/")
  } else {
    type <- paste0(type, "2-step floating catchment area")
  }
  wd <- crossprod(w, cv)
  wd[abs(wd) < .Machine$double.eps] <- 1
  r <- as.numeric(w %*% (pv / wd))
  if (!missing(return_type)) {
    return_type <- tolower(substr(return_type, 1, 1))
    if ("n" == return_type) {
      type <- paste(type, "(normalized)")
      r <- r / (sum(r * cv) / sum(cv))
    } else if ("r" == return_type) {
      type <- paste(type, "(resources per region)")
      r <- r * cv
    }
  } else {
    type <- paste(type, "(resources per consumer)")
  }
  if (print_type) cli_bullets(c(v = paste("calculated", type)))
  r
}
