test_that("example results align with manual", {
  demand <- c(5, 10, 50)
  supply <- c(5, 10)
  cost <- matrix(c(5, 35, 40, 55, 15, 25), ncol = 2)

  weight <- cost < 30
  access <- as.numeric(weight %*% (supply / crossprod(weight, demand)))
  expect_equal(catchment_ratio(demand, supply, cost, 30), access)

  # kernel density 2-step floating catchment area
  weight <- exp(-cost^2 / 8)
  access <- as.numeric(weight %*% (supply / crossprod(weight, demand)))
  expect_equal(catchment_ratio(demand, supply, cost, "gaussian"), access)
  expect_equal(catchment_ratio(demand, supply, cost, "dnorm"), access)
  expect_equal(catchment_ratio(demand, supply, cost, function(cost) dnorm(cost, 0, 2)), access)
  expect_equal(catchment_ratio(demand, supply, cost, weight), access)

  # enhanced 2-step floating catchment area
  weight <- (cost <= 60) * .22
  weight[cost <= 40] <- .68
  weight[cost <= 20] <- 1
  access <- as.numeric(weight %*% (supply / crossprod(weight, demand)))
  expect_equal(catchment_ratio(demand, supply, cost, list(c(60, .22), c(40, .68), c(20, 1))), access)

  # 3-step floating catchment area
  weight <- weight * weight / rowSums(weight)
  access <- as.numeric(weight %*% (supply / crossprod(weight, demand)))
  expect_equal(catchment_ratio(demand, supply, cost, list(c(60, .22), c(40, .68), c(20, 1)), TRUE), access)
})
