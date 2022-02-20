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
  expect_identical(connections, catchment_connections(
    sf::st_coordinates(pop$consumers),
    sf::st_coordinates(pop$providers)
  ))
  expect_identical(connections, catchment_connections(
    cbind(pop$consumers, sf::st_coordinates(pop$consumers)),
    cbind(pop$providers, sf::st_coordinates(pop$providers))
  ))
  pop$consumers$id <- seq_len(nrow(pop$consumers))
  expect_identical(connections, catchment_connections(pop$consumers, pop$providers, from_id = "id", to_id = "id"))
  c_rid <- sample(pop$consumers$id)
  p_rid <- sample(pop$providers$id)
  r_connections <- catchment_connections(
    pop$consumers[c_rid, ], pop$providers[p_rid, ],
    from_id = c_rid, to_id = p_rid
  )
  r_connections <- r_connections[order(as.numeric(r_connections$from), as.numeric(r_connections$to)), ]
  rownames(r_connections) <- NULL
  expect_identical(connections, r_connections)
  expect_identical(connections, catchment_connections(pop$consumers, pop$providers, cost = pop$cost))
  expect_identical(connections, catchment_connections(pop$consumers, pop$providers, cost = pop$cost[, p_rid]))
  rownames(pop$cost) <- pop$consumers$id
  expect_identical(connections, catchment_connections(pop$consumers, pop$providers, cost = pop$cost[c_rid, ]))
  expect_identical(connections, catchment_connections(pop$consumers, pop$providers, cost = as.numeric(pop$cost)))
  expect_identical(connections, catchment_connections(pop$consumers, pop$providers, cost = pop$cost, weight = pop$weight))
  expect_identical(connections, catchment_connections(pop$consumers, pop$providers, weight = pop$weight))
  expect_identical(connections, catchment_connections(pop$consumers, pop$providers, weight = as.numeric(pop$weight)))
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
