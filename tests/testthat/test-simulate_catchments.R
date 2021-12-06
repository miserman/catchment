test_that("returned results can be replicated", {
  pop <- simulate_catchments()
  expect_identical(pop$consumers$access, catchment_ratio(pop$consumers, pop$providers, pop$cost, pop$weight))
})
