test_that("connections line up as expected", {
  pop <- simulate_catchments()
  connections <- catchment_connections(pop$consumers, pop$providers)
  manual <- do.call(rbind, lapply(seq_len(nrow(pop$consumers)), function(r) {
    data.frame(from = r, to = which(pop$weight[r, ] != 0))
  }))
  expect_identical(as.integer(connections$from), manual$from)
  expect_identical(as.integer(connections$to), manual$to)
})
