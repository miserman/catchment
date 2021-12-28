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
#' Costs equal to \code{0} are treated as missing, so any truly \code{0} costs should be set to some minimal value.
#' @param weight Means of defining catchment areas and their topology / friction / impedance. The simplest is a single
#' number representing a maximum distance between \code{consumers} and \code{providers} (2-step floating catchment area;
#' 2SFCA; Luo & Wang, 2003).
#' An enhancement of this is a list of vectors with two values each: the first is a distance, and the second a weight
#' to associate with that distance (e.g., \code{list(c(10, 1), c(20, .5))}, which will give consumers within a
#' \code{cost} of 10 full weight, and those within a \code{cost} of 20 half weight; enhanced 2-step floating catchment
#' area; E2SFCA; Lou & Qi 2009). If a character, refers to a weighting function (kernel density 2-step floating
#' catchment area; KD2SFCA; Dai, 2010; in order from most gradual to steepest between costs of \code{1} and \code{6}):
#' \tabular{ll}{
#'   \code{linear} (\code{li}) \tab \code{{w <- (scale - cost) / scale; w[w < 0] <- 0; w}}\cr
#'   \code{gaussian} (\code{ga}) \tab \code{exp(-cost^2 / (2 * scale ^ 2))}\cr
#'   \code{d*} (name of a density function; e.g., \code{"dnorm"}) \tab
#'     \code{weight(cost, 0, scale)}
#'   \cr
#'   \code{p*} (name of a distribution function; e.g., \code{"pnorm"}) \tab \code{weight(-cost, 0, scale)}\cr
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
#' and 3-step floating catchment areas (3SFCA; Wan, Zou, & Sternberg, 2012).
#' @param scale Numeric scaling factor if \code{weight} is the name of a decay function.
#' @param max_cost Numeric limit on \code{cost}. This is the same as setting \code{weight} to a single value,
#' or specifying a list of steps as \code{weight} (where the most distant step is effectively \code{max_cost}),
#' although a single-value weight is exclusive (\code{cost < weight}) where steps are inclusive. This is most useful
#' when \code{weight} is a weighing function, where \code{max_cost} will trim the tail of the weight distribution.
#' @param adjust_consumers,adjust_providers A function to adjust weights when applied to \code{consumers} or
#' \code{providers}; should take the sparse weight matrix as its first argument, and return an adjusted matrix of the
#' same type. For example, you could square provider weights for the modified 2-step floating catchment area
#' (M2SFCA; Delamater, 2013) with \code{adjust_providers = function(w) w ^ 2}, or standardize both weights for the
#' balanced floating catchment area (BFCA; Paez, Higgins, & Vivona, 2019) with \code{adjust_consumers = }
#' \code{function(w) w / rowSums(w)} and \code{adjust_providers = function(w) sweep(w, 2, colSums(w), "/")}.
#' When weights are adjusted independently in this way, region scores will likely no longer sum to the sum
#' of \code{providers} (fewer than the total number of providers will be distributed).
#' @param return_type Determines the values that are returned: \code{"original"} (default) for \code{providers}
#' per \code{consumers} (e.g., how many, likely fractional, doctors are accessible by each person within each region),
#' \code{"region"} for number of \code{providers} per \code{consumers} entry (\code{consumers * original}; e.g.,
#' how many doctors are accessible within each region), or \code{"normalized"} for \code{original} divided by
#' \code{sum(region) / sum(consumers)}.
#' @param consumers_commutes A square, consumers source x consumers origin matrix with counts of origins,
#' used to specify multiple possible origins for each consumer location (e.g., consumers living in location 1
#' may work in locations 1 and 3, so the first row of \code{consumers_commutes} should have values in columns 1 and 3).
#' This can also be entered in place of \code{consumers}, assuming it includes all consumers (e.g., in a worker commute
#' matrix, you may need to add non-workers to the diagonal, if they are also consumers).
#' @param consumers_id,consumers_value,consumers_location,providers_id,providers_value,providers_location Column
#' names in \code{consumers} and/or \code{providers} to extract IDs, values, and location data (referring to a single
#' \code{sf} geometry column, or multiple columns with coordinates). These can also be used to directly enter
#' ID, value, and/or location vectors (or matrices for location coordinates).
#' @param verbose Logical; if \code{TRUE}, will print logs, and the type of floating catchment area that was calculated.
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
#' step_weights <- list(c(60, .22), c(40, .68), c(20, 1))
#' catchment_ratio(pop, doc, travel_time, step_weights)
#'
#' # modified 2-step floating catchment area
#' catchment_ratio(pop, doc, travel_time, step_weights, adjust_providers = function(m) m^2)
#'
#' # balanced 2-step floating catchment area
#' catchment_ratio(
#'   pop, doc, travel_time, step_weights,
#'   adjust_consumers = function(w) sweep(w, 1, rowSums(w), "/", FALSE),
#'   adjust_providers = function(w) sweep(w, 2, colSums(w), "/", FALSE),
#' )
#'
#' # 3-step floating catchment area
#' catchment_ratio(pop, doc, travel_time, step_weights, normalize_weight = TRUE)
#'
#' # visualized weight functions
#' if (require("splot", quietly = TRUE)) {
#'   cost <- 1:10
#'   scale <- 2
#'   splot(list(
#'     linear = (10 - cost) / 10,
#'     gaussian = exp(-cost^2 / (2 * scale^2)),
#'     dnorm = dnorm(cost, 0, scale),
#'     pnorm = pnorm(-cost, 0, scale),
#'     gravity = sqrt(1 / cost^scale),
#'     logarithmic = 1 / (1 + log(cost, scale)),
#'     logistic = 1 / (1 + exp(scale * cost)),
#'     exponential = exp(-cost * scale)
#'   ) ~ cost, title = "Decay Functions", laby = "Weight", labx = "Cost", lines = "con", note = FALSE)
#' }
#' @return \code{catchment_ratio}: A vector with an access score (determined by \code{return_type})
#' for each entry in \code{consumers}.
#' @references
#' Dai, D. (2010). Black residential segregation, disparities in spatial access to health care facilities, and
#' late-stage breast cancer diagnosis in metropolitan Detroit. \emph{Health & place, 16}, 1038-1052.
#' doi: \href{https://doi.org/10.1016/j.healthplace.2010.06.012}{10.1016/j.healthplace.2010.06.012}
#'
#' Delamater, P. L. (2013). Spatial accessibility in suboptimally configured health care systems: a modified
#' two-step floating catchment area (M2SFCA) metric. \emph{Health & place, 24}, 30-43.
#' doi: \href{https://doi.org/10.1016/j.healthplace.2013.07.012}{10.1016/j.healthplace.2013.07.012}
#'
#' Lou, W. & Qi, Y. (2009). An enhanced two-step floating catchment area (E2SFCA) method for measuring spatial
#' accessibility to primary care physicians. \emph{Health & Place, 15}, 1100-1107.
#' doi: \href{https://doi.org/10.1016/j.healthplace.2009.06.002}{10.1016/j.healthplace.2009.06.002}
#'
#' Luo, W. & Wang, F. (2003). Measures of spatial accessibility to health care in a GIS environment: synthesis
#' and a case study in the Chicago region. \emph{Environment and Planning B: Planning and Design, 30}, 865-884.
#' doi: \href{https://doi.org/10.1068/b29120}{10.1068/b29120}
#'
#' Paez, A., Higgins, C. D., & Vivona, S. F. (2019). Demand and level of service inflation in Floating Catchment
#' Area (FCA) methods. \emph{Plos one, 14}, e0218773. doi:
#' \href{https://doi.org/10.1371/journal.pone.0218773}{10.1371/journal.pone.0218773}
#'
#' Wan, N., Zou, B., & Sternberg, T. (2012). A three-step floating catchment area method for analyzing spatial
#' access to health services. \emph{International Journal of Geographical Information Science, 26}, 1073-1089.
#' doi: \href{https://doi.org/10.1080/13658816.2011.624987}{10.1080/13658816.2011.624987}
#' @export

catchment_ratio <- function(consumers = NULL, providers = NULL, cost = NULL, weight = NULL, normalize_weight = FALSE,
                            scale = 2, max_cost = NULL, adjust_consumers = NULL, adjust_providers = NULL,
                            return_type = "original", consumers_commutes = NULL, consumers_id = "GEOID",
                            consumers_value = "count", consumers_location = c("X", "Y"), providers_id = "GEOID",
                            providers_value = "count", providers_location = c("X", "Y"), verbose = FALSE) {
  if (verbose) cli_rule("Calculating a Floating Catchment Area Ratio")
  type <- ""
  if (is.null(consumers_commutes)) {
    if (!missing(consumers_value)) {
      dims <- dim(consumers_value)
      if (!is.null(dims) && dims[1] == dims[2]) {
        consumers_commutes <- consumers_value
        if (is.null(consumers)) consumers <- rowSums(consumers_commutes)
      }
    } else {
      dims <- dim(consumers)
      if (missing(consumers_id) && missing(consumers_location) && !is.null(dims) && dims[1] == dims[2]) {
        consumers_commutes <- consumers
        consumers <- rowSums(consumers)
      }
    }
  }
  input_data <- c(is.null(dim(providers)), is.null(dim(consumers)))
  # getting provider and consumer value vectors
  cv <- if (input_data[2]) {
    if (is.null(consumers)) {
      if (is.numeric(consumers_value)) {
        consumers_value
      } else {
        cli_abort("{.arg consumers} is not specified, and {.arg consumers_value} is non-numeric")
      }
    } else {
      if (verbose) cli_alert_info("consumers value: {.arg consumers} vector")
      as.numeric(consumers)
    }
  } else if (length(consumers_value) == 1 && consumers_value %in% colnames(consumers)) {
    if (verbose) cli_alert_info("consumers value: {.var {consumers_value}} column")
    as.numeric(consumers[, consumers_value, drop = TRUE])
  } else {
    if (verbose) cli_alert_info("consumers value: {.field 1}")
    rep(1, nrow(consumers))
  }
  if (!length(cv)) cli_abort("failed to recognize values in {.arg consumers}")
  pv <- if (input_data[1]) {
    if (is.null(providers)) {
      if (is.numeric(providers_value)) {
        providers_value
      } else {
        cli_abort("{.arg providers} is not specified, and {.arg providers_value} is non-numeric")
      }
    } else {
      if (verbose) cli_alert_info("providers value: {.arg providers} vector")
      as.numeric(providers)
    }
  } else if (length(providers_value) == 1 && providers_value %in% colnames(providers)) {
    if (verbose) cli_alert_info("providers value: {.var {providers_value}} column")
    as.numeric(providers[, providers_value, drop = TRUE])
  } else {
    if (verbose) cli_alert_info("providers value: {.field 1}")
    rep(1, nrow(providers))
  }
  if (!length(pv)) cli_abort("failed to recognize values in {.arg providers}")
  # getting provider and consumer ids
  cid <- if (input_data[1]) {
    if (!is.null(names(consumers))) {
      if (verbose) cli_alert_info("consumers id: {.arg consumers} names")
      names(consumers)
    } else if (!missing(consumers_id) && length(cv) == length(consumers_id)) {
      if (verbose) cli_alert_info("consumers id: {.arg consumers_id} vector")
      consumers_id
    } else if (!is.null(colnames(cost))) {
      if (verbose) cli_alert_info("consumers id: {.arg cost} column names")
      colnames(cost)
    } else {
      if (verbose) cli_alert_info("consumers id: sequence along consumers value")
      seq_along(cv)
    }
  } else if (length(consumers_id) == 1 && consumers_id %in% colnames(consumers)) {
    if (verbose) cli_alert_info("consumers id: {.var {consumers_id}} column")
    consumers[, consumers_id, drop = TRUE]
  } else {
    if (verbose) cli_alert_info("consumers id: sequence along consumers value")
    seq_along(cv)
  }
  if (!length(cid)) cli_abort("failed to recognize IDs in {.arg consumers}")
  pid <- if (input_data[1]) {
    if (!is.null(names(providers))) {
      if (verbose) cli_alert_info("providers id: {.arg providers} names")
      names(providers)
    } else if (!missing(providers_id) && length(pv) == length(providers_id)) {
      if (verbose) cli_alert_info("providers id: {.arg providers_id} vector")
      providers_id
    } else if (!is.null(colnames(cost))) {
      if (verbose) cli_alert_info("providers id: {.arg cost} column names")
      colnames(cost)
    } else {
      if (verbose) cli_alert_info("providers id: sequence along providers value")
      seq_along(pv)
    }
  } else if (length(providers_id) == 1 && providers_id %in% colnames(providers)) {
    if (verbose) cli_alert_info("providers id: {.var {providers_id}} column")
    providers[, providers_id, drop = TRUE]
  } else {
    if (verbose) cli_alert_info("providers id: sequence along providers value")
    seq_along(pv)
  }
  if (!length(pid)) cli_abort("failed to recognize IDs in {.arg providers}")
  if (is.null(cost)) {
    if (missing(consumers_location) && !all(consumers_location %in% colnames(consumers))) {
      consumers_location <- "geometry"
    }
    if (missing(providers_location) && !all(providers_location %in% colnames(providers))) {
      providers_location <- "geometry"
    }
    if ((!is.character(providers_location) || all(providers_location %in% colnames(providers))) &&
      (!is.character(consumers_location) || all(consumers_location %in% colnames(consumers)))) {
      if (verbose) cli_bullets("calculating cost from locations...")
      ccords <- if (is.character(consumers_location)) {
        if (verbose) cli_alert_info("consumers location: {.arg consumers} columns ({.field {consumers_location}})")
        consumers[, consumers_location, drop = TRUE]
      } else {
        if (verbose) cli_alert_info("consumers location: {.arg consumers_location}")
        consumers_location
      }
      if (any(grepl("^sf", class(ccords)))) {
        if (ncol(consumers_location) == 2 && is.numeric(consumers_location[, 1]) && is.numeric(consumers_location[, 2])) {
          if (verbose) cli_alert_info("dropping {.arg consumers_location} geometry")
          consumers_location <- consumers_location[, 1:2, drop = TRUE]
        } else {
          if (any(grepl("POLY", class(ccords), fixed = TRUE))) {
            if (verbose) cli_alert_info("calculating centroids of consumers location geometry")
            ccords <- st_centroid(ccords)
          }
          if (verbose) cli_alert_info("using coordinates from consumers location geometry")
          ccords <- st_coordinates(ccords)
        }
      }
      pcords <- if (is.character(providers_location)) {
        if (verbose) cli_alert_info("providers location: {.arg providers} columns ({.field {providers_location}})")
        providers[, providers_location, drop = FALSE]
      } else {
        if (verbose) cli_alert_info("providers location: {.arg providers_location}")
        providers_location
      }
      if (any(grepl("^sf", class(pcords)))) {
        if (ncol(providers_location) == 2 && is.numeric(providers_location[, 1]) && is.numeric(providers_location[, 2])) {
          if (verbose) cli_alert_info("dropping {.arg providers_location} geometry")
          providers_location <- providers_location[, 1:2, drop = FALSE]
        } else {
          if (any(grepl("POLY", class(pcords), fixed = TRUE))) {
            if (verbose) cli_alert_info("calculating centroids of providers location geometry")
            pcords <- st_centroid(pcords)
          }
          if (verbose) cli_alert_info("using coordinates from providers location geometry")
          pcords <- st_coordinates(pcords)
        }
      }
      if (ncol(pcords) == ncol(ccords)) {
        if (verbose) cli_alert_info("cost: calculated Euclidean distances")
        cost <- 1 / lma_simets(ccords, pcords, metric = "euclidean", pairwise = TRUE) - 1
        if (is.null(dim(cost))) {
          cost <- t(cost)
          if (nrow(pcords) == 1) cost <- t(cost)
        }
      } else {
        cli_abort(
          "{.arg cost} is NULL, and failed to calculate it from provided locations (differing number of columns)"
        )
      }
    }
    if (is.null(cost)) {
      if (verbose) cli_alert_info("cost: {.feild 1}")
      cost <- sparseMatrix(
        {},
        {},
        dims = c(length(cv), length(pv)),
        dimnames = list(cid, pid)
      )
      cost[] <- TRUE
    }
  } else if (is.null(dim(cost))) {
    cost <- t(cost)
    if (length(pv) == 1) cost <- t(cost)
    if (identical(dim(cost), c(length(cv), length(pv)))) {
      if (verbose) cli_alert_info("cost: {.arg cost} vector")
    } else {
      cli_abort("{.arg cost} must be a matrix-like object")
    }
  } else if (verbose) cli_alert_info("cost: {.arg cost} matrix")
  if (is.data.frame(cost)) cost <- as.matrix(cost)
  w <- catchment_weight(cost, weight, max_cost = max_cost, scale = scale, normalize_weight = FALSE, verbose = verbose)
  if (!is.null(consumers_commutes)) {
    dims <- dim(consumers_commutes)
    if (dims[1] != dims[2]) cli_abort("{.arg consumers_commutes} must be a square matrix")
    if (dims[1] != length(cv)) {
      if (all(cid %in% colnames(consumers_commutes)) || (is.numeric(cid) && dims[1] <= max(cid))) {
        consumers_commutes <- consumers_commutes[cid, cid]
      } else {
        cli_abort("{.arg consumers_commutes} could not be resolved with {.arg consumers}")
      }
    }
    noncommuters <- diag(consumers_commutes) / cv
    if (any(cv == 0)) noncommuters[cv == 0] <- 1
    diag(consumers_commutes) <- 0
    diag(consumers_commutes) <- rowSums(consumers_commutes)
    consumers_commutes <- consumers_commutes / rowSums(consumers_commutes)
    w_total <- rowSums(w)
    w_commutes <- consumers_commutes %*% (w / w_total) * w_total
    w <- w_commutes * (1 - noncommuters) + weight * noncommuters
  }
  wr <- rowSums(w)
  if (!is.null(rownames(w)) && !all(cid == rownames(w)) && all(cid %in% rownames(w))) {
    if (verbose) cli_alert_info("selected weight rows by consumers id names")
    w <- w[cid, ]
  }
  if (!is.null(colnames(w)) && !all(pid == colnames(w)) && all(pid %in% colnames(w))) {
    if (verbose) cli_alert_info("selected weight columns by providers id names")
    w <- w[, pid]
  }
  if (nrow(w) != length(cv) && is.numeric(cid) && min(cid) > 0 && max(cid) <= nrow(w)) {
    if (verbose) cli_alert_info("selected weight rows by consumers id indices")
    w <- w[cid, ]
  }
  if (ncol(w) != length(pv) && is.numeric(pid) && min(pid) > 0 && max(pid) <= ncol(w)) {
    if (verbose) cli_alert_info("selected weight rows by providers id indices")
    w <- w[, pid]
  }
  if (!all(dim(w) == c(length(cv), length(pv)))) {
    cli_abort("failed to align weight matrix with consumer and provider values")
  }
  if (normalize_weight) {
    # adjust weights by selection probability (3-step){
    if (verbose) cli_alert_info("normalizing weight")
    type <- paste0(type, "3-step floating catchment area")
    if (!is.null(rownames(w)) && all(rownames(w) %in% names(wr))) wr <- wr[rownames(w)]
    wr[wr == 0] <- 1
    w <- w * sweep(w, 1, wr, "/", FALSE)
  } else {
    type <- paste0(type, "2-step floating catchment area")
  }
  wd <- crossprod(if (is.function(adjust_consumers)) {
    environment(adjust_consumers) <- environment()
    adjust_consumers(w)
  } else {
    w
  }, cv)
  wd[abs(wd) < .Machine$double.eps] <- 1
  r <- as.numeric((if (is.function(adjust_providers)) {
    environment(adjust_providers) <- environment()
    adjust_providers(w)
  } else {
    w
  }) %*% (pv / wd))
  if (!missing(return_type)) {
    return_type <- tolower(substr(return_type, 1, 1))
    if ("n" == return_type) {
      type <- paste(type, "(normalized)")
      if (verbose) cli_alert_info("normalizing ratios")
      r <- r / (sum(r * cv) / sum(cv))
    } else if ("r" == return_type) {
      type <- paste(type, "(resources per region)")
      if (verbose) cli_alert_info("multiplying ratios by consumers value")
      r <- r * cv
    }
  } else {
    type <- paste(type, "(resources per consumer)")
  }
  if (verbose) cli_bullets(c(v = paste("calculated", type)))
  r
}
