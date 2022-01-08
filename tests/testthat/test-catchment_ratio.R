test_that("example results align with manual", {
  demand <- c(5, 10, 50)
  supply <- c(5, 10)
  cost <- matrix(c(5, 35, 40, 55, 15, 25), ncol = 2)

  weight <- cost < 30
  access <- as.numeric(weight %*% (supply / crossprod(weight, demand)))
  expect_equal(catchment_ratio(demand, supply, cost, 30), access)
  expect_equal(catchment_ratio(demand, supply, cost, 30, return_type = "region"), access * demand)
  expect_equal(
    catchment_ratio(demand, supply, cost, 30, return_type = "norm"),
    access / (sum(access * demand) / sum(demand))
  )

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

test_that("alternative inputs work", {
  library(lingmatch)
  demand <- data.frame(
    id = as.character(1:50),
    d = sample(10:50, 50, TRUE),
    lat = rnorm(50),
    long = rnorm(50)
  )
  rownames(demand) <- demand$id
  supply <- data.frame(
    id = as.character(100:109),
    s = sample(0:10, 10, TRUE),
    lat = rnorm(10),
    long = rnorm(10)
  )
  rownames(supply) <- supply$id
  cost <- 1 / lma_simets(demand[, c("lat", "long")], supply[, c("lat", "long")], "euc", pairwise = TRUE) - 1
  m <- catchment_ratio(
    demand, supply, cost,
    weight = 1,
    consumers_id = "id", consumers_value = "d",
    providers_id = "id", providers_value = "s"
  )
  expect_equal(m, catchment_ratio(
    demand, supply[sample(supply$id), ], cost,
    weight = 1,
    consumers_id = "id", consumers_value = "d",
    providers_id = "id", providers_value = "s"
  ))
  expect_equal(m, catchment_ratio(
    as.matrix(demand), as.matrix(supply), as.data.frame(as.matrix(cost)),
    weight = 1,
    consumers_id = "id", consumers_value = "d",
    providers_id = "id", providers_value = "s"
  ))
  expect_equal(m, catchment_ratio(
    demand, supply,
    weight = 1,
    consumers_id = "id", consumers_value = "d", consumers_location = c("lat", "long"),
    providers_id = "id", providers_value = "s", providers_location = c("lat", "long")
  ))
  expect_equal(m, catchment_ratio(
    cost = cost, weight = 1,
    consumers_id = demand$id, consumers_value = demand$d,
    providers_id = supply$id, providers_value = supply$s
  ))
  expect_equal(m, catchment_ratio(
    cost = cost, weight = 1,
    consumers_id = demand$id, consumers_value = demand$d,
    providers_id = supply$id, providers_value = supply$s
  ))
  expect_equal(m, catchment_ratio(
    weight = 1,
    consumers_location = demand[, c("lat", "long")],
    providers_location = supply[, c("lat", "long")],
    consumers_id = demand$id, consumers_value = demand$d,
    providers_id = supply$id, providers_value = supply$s
  ))
})

test_that("consumers and providers get aligned with cost and weights", {
  co <- structure(sample(100:5000, 100), names = paste0("c", 1:100))
  pr <- structure(sample(1:10, 50, TRUE), names = paste0("p", 1:50))
  cost <- matrix(rpois(5e3, 4) * 10, 100, 50, dimnames = list(names(co), names(pr)))

  ex <- 5:7
  cost[, ex] <- 99
  full <- catchment_ratio(co, pr, cost, 40)
  expect_equal(full, catchment_ratio(co, pr, cost[, -ex], 40))
  expect_equal(full, catchment_ratio(co, pr[-ex], cost, 40))
  cost[ex, ] <- 0
  full <- catchment_ratio(co, pr, cost, 40)
  expect_equal(full, catchment_ratio(co, pr, cost[-ex, ], 40))
  expect_equal(full, catchment_ratio(co, pr[-ex], cost, 40))

  is <- c(1, 5, 10)
  full <- catchment_ratio(co, pr[is], cost[, is], 40)
  expect_equal(full, catchment_ratio(co, pr, cost[, is], 40))
  expect_equal(full, catchment_ratio(co, pr[is], cost, 40))
  full <- catchment_ratio(co[is], pr, cost[is, ], 40)
  expect_equal(full, catchment_ratio(co, pr, cost[is, ], 40)[is])
  expect_equal(full, catchment_ratio(co[is], pr, cost, 40))
})
