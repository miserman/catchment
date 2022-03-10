#' Aggregate Floating Catchment Area Ratios
#'
#' Calculate aggregate sums or means, potentially weighted by number of consumers.
#'
#' @param from A matrix-like object with IDs, values, and number of consumers, or a vector
#' of values. These are the lower-level entities to be aggregated up to entries in \code{to}.
#' @param to A matrix-like object with IDs and total consumers. These are higher-level
#' entities containing \code{from} entities.
#' @param id The column name of IDs in \code{from}, or a vector of IDs.
#' @param value The column name of values in \code{from} or a vector of values to be aggregated.
#' @param consumers The column name of consumer totals in \code{from}, or a vector
#' of said totals of the same length as \code{from} IDs.
#' @param to_id The column name of IDs in \code{to}, or a vector if IDs.
#' @param to_consumers The column name of consumer totals in \code{to}, or a vector
#' of said totals of the same length as \code{to} IDs.
#' @param map A named list, with names as or corresponding to \code{to} IDs, and vectors
#' containing \code{from} IDs as entries.
#' @param original_from Logical indicating whether the \code{from} values are original (default).
#' If \code{from} values are original and \code{consumers} are specified, \code{from} values will be
#' multiplied by \code{consumers} before aggregation. Set this to \code{FALSE} if this has already been done
#' (such as if you set \code{return_type} to \code{"region"} in \code{\link{catchment_ratio}}).
#' @param return_type Determines the values that are returned: \code{"original"} (default) for an average
#' for each \code{to} ID based on the region sum of \code{consumers}, or specified \code{to_consumers}.
#' If a number, original scores will be multiplied by the number then averaged like \code{"original"}. If
#' \code{"region"}, the sum for each \code{to} ID is returned.
#' @param verbose Logical; if \code{TRUE}, will print a log of the aggregation process.
#' @examples
#' # lower-level entries prepended with higher-level IDs (like GEOIDs)
#' higher <- data.frame(id = c(100, 110), pop = c(500, 500))
#' lower <- data.frame(
#'   id = paste0(higher$id, 1:5),
#'   value = 1:5 / 10,
#'   pop = c(50, 10, 100, 40, 30)
#' )
#' catchment_aggregate(
#'   lower, higher,
#'   consumers = "pop", id = "id", value = "value",
#'   to_id = "id", verbose = TRUE
#' )
#'
#' # same thing with IDs specified in a map
#' catchment_aggregate(
#'   lower, higher,
#'   consumers = "pop", id = "id", value = "value",
#'   map = list("100" = c(1001, 1003, 1005), "110" = c(1102, 1104)), verbose = TRUE
#' )
#' @return A vector with an aggregate value (determined by \code{return_type})
#' for each entry in \code{to}.
#' @export

catchment_aggregate <- function(from, to = NULL, id = "GEOID", value = "access", consumers = "population",
                                to_id = id, to_consumers = consumers, map = NULL, original_from = TRUE,
                                return_type = "original", verbose = FALSE) {
  if (verbose) cli_rule("Aggregating Values")
  if (missing(from)) {
    if (!missing(value)) {
      from <- value
    } else {
      cli_abort("{.arg from} must be provided")
    }
  }
  fid <- if (length(id) > 1) {
    if (verbose) cli_alert_info("from IDs: {.arg id} vector")
    id
  } else if (is.null(dim(from))) {
    if (is.null(names(from))) cli_abort("{.arg id} must be specified if they are not provided in {.arg from}")
    if (verbose) cli_alert_info("from IDs: {.arg from} names")
    names(from)
  } else {
    if (verbose) cli_alert_info("from IDs: {.var {id}} column")
    from[[id]]
  }
  if (is.null(fid)) cli_abort("failed to resolve from IDs")
  fv <- if (length(value) > 1) {
    if (verbose) cli_alert_info("from values: {.arg value} vector")
    value
  } else if (is.null(dim(from))) {
    if (verbose) cli_alert_info("from values: entered vector")
    from
  } else if (!is.null(value) && value %in% colnames(from)) {
    if (verbose) cli_alert_info("from values: {.var {value}} column")
    from[[value]]
  } else {
    if (verbose) cli_alert_info("from values: {.feild 1}")
    rep(1, length(fid))
  }
  if (is.null(fv) || length(fv) != length(fid)) cli_abort("failed to resolve {.arg from} values")

  n <- length(fv)
  if (length(fid) != n) cli_abort("IDs and values were not the same length")
  cv <- if (length(consumers) > 1) {
    if (verbose) cli_alert_info("consumers: {.arg consumers} vector")
    consumers
  } else if (!is.null(consumers) && consumers %in% colnames(from)) {
    if (verbose) cli_alert_info("consumers: {.var {consumers}} column")
    from[[consumers]]
  } else {
    NULL
  }
  if (!is.null(cv) && length(cv) != n) cli_abort("{.arg consumers} was not the same length as values")

  if (original_from && !is.null(cv) && length(cv) == n) fv <- fv * cv
  # `from` values should now be either per-region values or normalized

  # mapping `from` ids to `to` ids
  if (is.null(to)) to <- if (!missing(to_id) && (!is.character(to_id) || length(to_id) > 1)) to_id else from
  tid <- if (length(to_id) > 1) {
    if (verbose) cli_alert_info("to IDs: {.arg to_id} vector")
    to_id
  } else if (!is.null(to_id) && to_id %in% colnames(to)) {
    if (verbose) cli_alert_info("to IDs: {.var {to_id}} column")
    to[[to_id]]
  } else if (length(to) == 1) {
    if (is.function(to)) {
      if (verbose) cli_alert_info("to IDs: {.arg from} IDs converted by {.arg to} function")
      to(fid)
    } else if (is.numeric(to)) {
      if (verbose) cli_alert_info(paste("to IDs: first", to, "characters of {.arg from} IDs"))
      substr(fid, 1, to)
    } else {
      if (verbose) cli_alert_info("to IDs: {.var {to}} column of {.arg from}")
      from[[to]]
    }
  } else if (is.null(dim(to))) {
    if (verbose) cli_alert_info("to IDs: {.arg to} vector")
    to
  } else {
    NULL
  }
  if (is.null(tid)) cli_abort("failed to resolve {.arg to} IDs")

  tcv <- if (length(to_consumers) > 1) {
    if (verbose) cli_alert_info("to consumers: {.arg to_consumers} vector")
    to_consumers
  } else if (!is.null(to_consumers) && to_consumers %in% colnames(to)) {
    if (verbose) cli_alert_info("to consumers: {.var {to_consumers}} column")
    to[[to_consumers]]
  } else {
    NULL
  }
  aggregate_consumers <- is.null(tcv)
  fv <- cbind(fv, if (!aggregate_consumers || is.null(cv)) 1 else cv)

  if (is.list(map)) {
    if (is.null(names(map))) cli_abort("{.arg map} must have names")
    if (verbose) cli_alert_info("mapping: by map list")
    if (is.null(tid)) tid <- names(map)
    utid <- as.character(unique(tid))
    av <- vapply(utid, function(h) {
      l <- map[[h]]
      su <- fid %in% l
      if (is.null(l) || !length(su)) {
        return(c(NA, NA))
      }
      colSums(fv[su, , drop = FALSE], na.rm = TRUE)
    }, numeric(2))
  } else {
    if (length(tid) == n) {
      if (verbose) cli_alert_info("mapping: by {.arg to} ID")
      fid <- tid
      tid <- unique(tid)
    } else {
      tnc <- nchar(tid[1])
      if (tnc > nchar(fid[1]) || any(nchar(tid) != tnc)) {
        if (!any(fid %in% tid)) cli_abort("failed to figure out how to map {.arg from} ids to {.arg to} ids")
        if (verbose) cli_alert_info("mapping: by ID match")
      } else {
        if (verbose) cli_alert_info(paste("mapping: by first", tnc, "character match"))
        fid <- substr(fid, 1, tnc)
      }
    }
    utid <- as.character(unique(tid))
    av <- vapply(utid, function(h) {
      su <- fid == h
      if (any(su)) colSums(fv[su, , drop = FALSE], na.rm = TRUE) else c(NA, NA)
    }, numeric(2))
  }
  if (length(utid) != length(tid)) {
    if (verbose) cli_alert_info("to repeated IDs")
    av <- rbind(unname(av[1, ][tid]), unname(av[2, ][tid]))
  }

  if (is.numeric(return_type)) {
    if (verbose) cli_alert_info(paste("multiplying by", return_type))
    av[1, ] <- av[1, ] * return_type
  }
  if (is.character(return_type) && grepl("^r", return_type, TRUE)) {
    if (verbose) cli_bullets(c(v = "returning sum of value per {.arg to} ID"))
    av[1, ]
  } else {
    if (aggregate_consumers && !is.null(tcv) && ncol(av) != length(tcv)) aggregate_consumers <- FALSE
    if (verbose) {
      cli_bullets(c(v = paste0(
        "returning sum of value over total {.arg ",
        if (aggregate_consumers) "from" else "to",
        "} consumers per {.arg to} ID"
      )))
    }
    if (aggregate_consumers) {
      av[2, av[2, ] == 0] <- 1
      av[1, ] / av[2, ]
    } else {
      tcv[tcv == 0] <- 1
      av[1, ] / tcv
    }
  }
}
