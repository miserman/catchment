#' Extract Connections Within Catchment Areas
#'
#' Extract connected consumer (\code{from}) and provider (\code{to}) locations within a catchment area,
#' as defined by cost and weights.
#'
#' @param from Consumer locations; a matrix-like object with columns containing latitudes and longitudes
#' (identified by \code{from_coords}), and a column of IDs corresponding to rows in \code{cost} (identified by
#' \code{from_id}). Each identifying column can alternatively contain the associated matrix or vector, and
#' IDs can be in row names.
#' @param to Provider locations;a matrix-like object with columns containing latitudes and longitudes
#' (identified by \code{to_coords}), and a column of IDs corresponding to columns in \code{cost} (identified by
#' \code{to_id}). Each identifying column can alternatively contain the associated matrix or vector, and
#' IDs can be in row names.
#' @param cost A cost matrix, with row names corresponding to IDs of \code{from}, and column names corresponding to
#' IDs of \code{to}.
#' @param weight A weight matrix with the same dimensions as \code{cost}, such as returned from
#' \code{\link{catchment_weight}}.
#' @param ... Passes arguments to \code{catchment_weight} if \code{weight} is not a weight matrix.
#' @param return_type Specify whether to return a \code{data.frame} (default) with connection ids, weights, and costs,
#' an \code{sf} \code{data.frame} (\code{"sf"}) with linestring geometries added for each connection,
#' or a GeoJSON-formatted list (\code{"list"}).
#' @param from_id,from_coords,to_id,to_coords Names of ID and coordinate columns in \code{from} and \code{to},
#' or vectors of IDs and matrices of coordinates.
#' @examples
#' pop <- simulate_catchments()
#' connections <- catchment_connections(
#'   pop$consumers, pop$providers,
#'   weight = "gaussian", max_cost = 1,
#'   return_type = "sf"
#' )
#' if (require("leaflet", quiet = TRUE)) {
#'   leaflet() |>
#'     addPolylines(
#'       data = connections, weight = 3, color = "#777",
#'       highlightOptions = highlightOptions(color = "#fff", weight = 4),
#'       label = ~ paste0("From: ", from, ", To: ", to, ", Weight: ", weight, ", Cost: ", cost)
#'     ) |>
#'     addCircles(
#'       data = pop$consumers, label = ~ paste0("Consumers: ", count, ", Access: ", access)
#'     ) |>
#'     addCircles(
#'       data = pop$providers, color = "#000", label = ~ paste("Provider", id)
#'     )
#' }
#' @return An object with connection, weight, and cost information, with a format depending on \code{return_type}.
#' @export

catchment_connections <- function(from, to, cost = NULL, weight = 1, ..., return_type = "data.frame",
                                  from_id = "GEOID", from_coords = c("X", "Y"), to_id = "GEOID",
                                  to_coords = c("X", "Y")) {
  if (missing(from)) from <- from_coords
  if (missing(to)) to <- from_coords
  from_ids <- if (length(from_id) != 1) {
    from_id
  } else if (from_id %in% colnames(from)) {
    from[, from_id, drop = TRUE]
  } else {
    rownames(from)
  }
  if (is.null(from_ids)) cli_abort("failed to resolve {.arg from} ids")
  if (is.numeric(from_ids) && max(from_ids) > nrow(from)) from_ids <- as.character(from_ids)
  to_ids <- if (length(to_id) != 1) {
    to_id
  } else if (to_id %in% colnames(to)) {
    to[, to_id, drop = TRUE]
  } else {
    rownames(to)
  }
  if (is.null(to_id)) cli_abort("failed to resolve {.arg to} ids")
  if (is.numeric(to_ids) && max(to_ids) > nrow(to)) to_ids <- as.character(to_ids)
  if (is.character(from_coords) && all(from_coords %in% colnames(from))) {
    from <- from[, from_coords, drop = TRUE]
  } else if (any(grepl("^sfc?$", class(from)))) {
    from <- st_coordinates(from)
  } else if (!is.null(dim(from)) && ncol(from) > 1) {
    from <- from[, 1:2, drop = TRUE]
  }
  if (is.null(dim(from)) || ncol(from) != 2) cli_abort("{.arg from} must have 2 columns")
  if (!is.numeric(from[, 1]) || !is.numeric(from[, 2])) {
    cli_abort("{.arg from} columns must all be numeric")
  }
  if (length(from_ids) != nrow(from)) cli_abort("From IDs and coordinates do not align")
  if (is.character(to_coords) && all(to_coords %in% colnames(to))) {
    to <- to[, to_coords, drop = TRUE]
  } else if (any(grepl("^sfc?$", class(to)))) {
    to <- st_coordinates(to)
  } else if (!is.null(dim(to)) && ncol(to) > 1) {
    to <- to[, 1:2, drop = TRUE]
  }
  if (is.null(dim(to)) || ncol(to) != 2) cli_abort("{.arg to} must have 2 columns")
  if (!is.numeric(to[, 1]) || !is.numeric(to[, 2])) {
    cli_abort("{.arg to} columns must all be numeric")
  }
  if (length(to_ids) != nrow(to)) cli_abort("To IDs and coordinates do not align")
  if (is.null(cost)) cost <- 1 / lma_simets(from, to, "euclidean", symmetrical = TRUE) - 1
  if (is.null(dim(weight)) || length(substitute(...()))) {
    weight <- catchment_weight(cost = cost, weight = weight, ...)
  } else if (any(dim(cost) != dim(weight))) cli_abort("{.arg weight} does not align with {.arg cost}")
  fcoords <- split(from, from_ids)[from_ids]
  tcoords <- split(to, to_ids)[to_ids]
  if (grepl("^[slg]", return_type, TRUE)) {
    res <- list(
      type = "FeatureCollection",
      features = unlist(lapply(unname(which(rowSums(weight != 0) != 0)), function(r) {
        su <- which(weight[r, ] != 0)
        np <- names(su)
        w <- weight[r, su]
        wc <- cost[r, su]
        fc <- as.numeric(fcoords[[r]])
        tc <- tcoords[to_ids %in% np]
        lapply(seq_along(np), function(i) {
          list(
            type = "Feature",
            properties = list(from = from_ids[[r]], to = np[[i]], weight = w[[i]], cost = wc[[i]]),
            geometry = list(
              type = "LineString",
              coordinates = list(fc, as.numeric(tc[[np[[i]]]]))
            )
          )
        })
      }), FALSE)
    )
    if (!length(res$features)) cli_abort("there are no connections")
    if (grepl("^s", return_type, TRUE)) {
      read_sf(toJSON(res, auto_unbox = TRUE), as_tibble = FALSE)
    } else {
      res
    }
  } else {
    rs <- unname(which(rowSums(weight != 0) != 0))
    if (!length(rs)) cli_abort("there are no connections")
    do.call(rbind, lapply(rs, function(r) {
      su <- which(weight[r, ] != 0)
      data.frame(
        from = from_ids[[r]], to = names(su), weight = weight[r, su], cost = cost[r, su]
      )
    }))
  }
}
