#' Extract a network of consumers and providers
#'
#' Extract an interconnected set of consumers and providers from a set of connections within catchment areas.
#'
#' @param connections A list- or matrix-like object with \code{"from"} and \code{"to"} entries, such as that
#' returned from \code{\link{catchment_connections}}.
#' @param from_start,to_start The ID of a \code{from} (consumer) or \code{to} (provider) in \code{connections} from
#' which to trace out a network.
#' @examples
#' pop <- simulate_catchments()
#' connections <- catchment_connections(
#'   pop$consumers, pop$providers,
#'   weight = 1, max_cost = .5,
#' )
#' catchment_network(connections, 1)
#' @return A subsetted version of \code{connections} containing only connections traceable from \code{from_start}
#' and/or \code{to_start}.
#' @export

catchment_network <- function(connections, from_start = NULL, to_start = NULL) {
  if (any(!c("from", "to") %in% colnames(connections))) cli_abort("{.arg connections} not recognized")
  from <- connections[["from"]]
  to <- connections[["to"]]
  from_ids <- unique(from)
  froms <- if (is.null(from_start) && is.null(to_start)) from[1] else from_start
  tos <- to_start
  for (i in seq_along(from_ids)) {
    tos <- unique(c(tos, to[from %in% froms]))
    froms <- unique(c(froms, from[to %in% tos]))
    if (all(from_ids %in% froms)) break
  }
  connections[from %in% froms | to %in% tos, ]
}
