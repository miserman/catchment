test_that("all returned cases are connected", {
  pop <- simulate_catchments()
  connections <- catchment_connections(pop$consumers, pop$providers)

  network <- catchment_network(connections, 1)
  network <- table(network$from, network$to)
  expect_true(all(rowSums(network) != 0) && all(colSums(network) != 0))

  network <- catchment_network(connections, to_start = connections$to[1])
  network <- table(network$from, network$to)
  expect_true(all(rowSums(network) != 0) && all(colSums(network) != 0))
})
