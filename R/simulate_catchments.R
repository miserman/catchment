#' Simulate Catchment Areas
#'
#' Generate a random population of consumers and providers.
#'
#' @param n Number of consumers to generate.
#' @param weight Weight method used on Euclidean distances between consumers and providers.
#' @param ... Passes additional arguments to \code{\link{catchment_weight}}.
#' @param consumers_range A numeric vector with entries for minimal and maximal consumer counts
#' from which to sample.
#' @param long_center Longitude to randomly place consumers around.
#' @param long_spread Standard deviation of the consumer's longitude distribution.
#' @param lat_center Latitude to randomly place consumers around.
#' @param lat_spread Standard deviation of the consumer's latitude distribution.
#' @param as_sf Logical; if \code{FALSE}, consumer and provider data are returned in regular \code{data.frame}s.
#' @examples
#' pop <- simulate_catchments()
#' if (require("leaflet", quiet = TRUE)) {
#'   leaflet(pop$providers) |>
#'     addCircles(
#'       radius = 1e5, color = "#666", stroke = FALSE
#'     ) |>
#'     addCircles(
#'       data = pop$consumers, label = ~ paste0("Consumers: ", count, ", Access: ", access)
#'     ) |>
#'     addCircles(
#'       color = "#000", label = ~ paste("Provider", id)
#'     )
#' }
#' @return A list with \code{consumers}, \code{providers}, \code{cost}, and \code{weight} entries.
#' @export

simulate_catchments <- function(n = 100, weight = 1, ..., consumers_range = c(10, 5000), long_center = -99,
                                long_spread = 1, lat_center = 40, lat_spread = 1, as_sf = TRUE) {
  # roll consumers
  consumers <- data.frame(
    count = sample(seq(consumers_range[1], consumers_range[2]), n, TRUE),
    X = rnorm(n, long_center, long_spread),
    Y = rnorm(n, lat_center, lat_spread)
  )
  locs <- consumers[, c("X", "Y")]
  consumers_dist <- 1 / lma_simets(locs, "euclidean", symmetrical = TRUE) - 1
  diag(consumers_dist) <- Inf
  start <- which.max(diag(consumers_dist[, max.col(-consumers_dist)]))
  diag(consumers_dist) <- 0
  providers <- data.frame(id = seq_len(nrow(consumers)), X = 0, Y = 0)
  # placing providers
  for (i in providers$id) {
    sel <- c(start, which(consumers_dist[start, ] * catchment_weight(
      consumers_dist[start, , drop = FALSE],
      weight, ...
    ) != 0))
    providers[i, c("X", "Y")] <- colMeans(locs[sel, , drop = FALSE] + rnorm(length(sel) * 2, 0, .1))
    r <- catchment_ratio(consumers, providers[seq_len(i), , drop = FALSE], weight = weight, ...)
    if (all(r != 0)) {
      providers <- providers[seq_len(i), ]
      break
    }
    start <- which.min(r)
  }
  consumers$access <- r
  cost <- 1 / lma_simets(consumers[, -1], providers[, -1], "euclidean", symmetrical = TRUE) - 1
  list(
    consumers = if (as_sf) st_as_sf(consumers, coords = c("X", "Y"), crs = 4326) else consumers,
    providers = if (as_sf) st_as_sf(providers, coords = c("X", "Y"), crs = 4326) else providers,
    cost = cost, weight = catchment_weight(cost, weight, ...)
  )
}
