test_that("connections line up as expected", {
  pop <- simulate_catchments()
  connections <- catchment_connections(pop$consumers, pop$providers)
  manual <- do.call(rbind, lapply(seq_len(nrow(pop$consumers)), function(r) {
    data.frame(from = r, to = which(pop$weight[r, ] != 0))
  }))
  expect_identical(as.integer(connections$from), manual$from)
  expect_identical(as.integer(connections$to), manual$to)
})

test_that("different inputs and outputs work", {
  pop <- simulate_catchments()
  connections <- catchment_connections(pop$consumers, pop$providers)
  expect_identical(connections, catchment_connections(from_coords = pop$consumers, to_coords = pop$providers))
  connections$cost <- round(connections$cost, 4)
  connections$weight <- as.integer(connections$weight)
  connnections_sf <- catchment_connections(pop$consumers, pop$providers, return_type = "sf")
  expect_identical(connections, connnections_sf[, 1:4, drop = TRUE])
  expect_identical(
    connnections_sf,
    sf::st_read(jsonlite::toJSON(
      catchment_connections(pop$consumers, pop$providers, return_type = "geojson"),
      auto_unbox = TRUE
    ), as_tibble = FALSE, quiet = TRUE)
  )
})
