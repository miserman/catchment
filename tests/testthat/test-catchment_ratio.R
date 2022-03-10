library(lingmatch)
library(sf)

test_that("example results align with manual", {
  demand <- c(5, 10, 50)
  supply <- c(5, 10)
  cost <- matrix(c(5, 35, 40, 55, 15, 25), ncol = 2)

  weight <- cost < 30
  access <- as.numeric(weight %*% (supply / crossprod(weight, demand)))
  expect_equal(catchment_ratio(demand, supply, cost, 30), access)
  expect_equal(catchment_ratio(demand, supply, cost, 30, return_type = "region"), access * demand)
  expect_equal(catchment_ratio(demand, supply, cost, 30, return_type = 1000), access * 1000)
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

  # commuter-based 2-step floating catchment area
  consumers_work <- matrix(abs(rnorm(length(demand) * length(demand))), length(demand))
  consumers_work <- consumers_work / rowSums(consumers_work) * demand
  nonworker_prop <- diag(consumers_work) / demand
  commuters <- consumers_work
  diag(commuters) <- 0
  diag(commuters) <- rowSums(commuters)
  commuters <- commuters / rowSums(commuters)
  consumer_weight <- rowSums(weight)
  weight_prop <- weight / consumer_weight
  weight_commute <- commuters %*% weight_prop * consumer_weight
  weight_commute <- weight_commute * (1 - nonworker_prop) + weight * nonworker_prop
  access <- catchment_ratio(demand, supply, weight = weight_commute)
  expect_equal(access, catchment_ratio(consumers_work, supply, weight = weight), tolerance = 1e13)
  expect_equal(access, catchment_ratio(consumers_value = consumers_work, providers_value = supply, weight = weight), tolerance = 1e13)
  expect_identical(access, catchment_ratio(demand, supply, weight = weight, consumers_commutes = consumers_work))

  names(demand) <- 1:3
  names(supply) <- 1:2
  dimnames(weight) <- list(1:3, 1:2)
  dimnames(consumers_work) <- list(1:3, 1:3)
  expect_identical(access, catchment_ratio(
    demand, supply,
    weight = weight,
    consumers_commutes = consumers_work[sample(1:3), sample(1:3)]
  ))

  # balanced 2-step floating catchment area
  access <- as.numeric(sweep(weight, 2, colSums(weight), "/", FALSE) %*% (supply / crossprod(sweep(weight, 1, rowSums(weight), "/", FALSE), demand)))
  expect_equal(
    access,
    catchment_ratio(
      demand, supply, cost, weight,
      adjust_consumers = function(w) sweep(w, 1, rowSums(w), "/", FALSE),
      adjust_providers = function(w) sweep(w, 2, colSums(w), "/", FALSE),
    )
  )

  # 3-step floating catchment area
  weight <- weight * weight / rowSums(weight)
  access <- as.numeric(weight %*% (supply / crossprod(weight, demand)))
  expect_equal(catchment_ratio(demand, supply, cost, list(c(60, .22), c(40, .68), c(20, 1)), TRUE), access)
})

test_that("alternative inputs work", {
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
    demand, supply, as.numeric(cost),
    weight = 1,
    consumers_id = "id", consumers_value = "d",
    providers_id = "id", providers_value = "s"
  ))
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

  demand <- st_as_sf(demand, coords = c("lat", "long"), remove = FALSE)
  supply <- st_as_sf(supply, coords = c("lat", "long"), remove = FALSE)
  expect_equal(m, catchment_ratio(
    demand, supply,
    weight = 1,
    consumers_id = "id", consumers_value = "d",
    providers_id = "id", providers_value = "s"
  ))
  expect_equal(m, catchment_ratio(
    demand, supply,
    weight = 1,
    consumers_id = "id", consumers_value = "d",
    providers_id = "id", providers_value = "s",
    consumers_location = c("lat", "long"),
    providers_location = c("lat", "long")
  ))
  expect_equal(m, catchment_ratio(
    demand, supply,
    weight = 1,
    consumers_id = "id", consumers_value = "d",
    providers_id = "id", providers_value = "s",
    consumers_location = demand[, c("lat", "long")],
    providers_location = supply[, c("lat", "long")]
  ))
})

test_that("consumers and providers get aligned with cost and weights", {
  co <- structure(sample(100:5000, 100), names = paste0("c", 1:100))
  pr <- structure(sample(1:10, 50, TRUE), names = paste0("p", 1:50))
  cost <- matrix(rpois(5e3, 4) * 10, 100, 50)

  ex <- 5:7
  cost[, ex] <- 99
  full <- catchment_ratio(co, pr, cost, 40)
  expect_equal(full, catchment_ratio(co, pr, cbind(cost, cost), 40, providers_id = 1:50))
  expect_equal(full, catchment_ratio(co, pr, rbind(cost, cost), 40, consumers_id = 1:100))

  dimnames(cost) <- list(1:100, 1:50)
  expect_equal(full, catchment_ratio(co, pr, cost[, -ex], 40, providers_id = 1:50))
  expect_equal(full, catchment_ratio(co, pr[-ex], cost, 40, providers_id = (1:50)[-ex], consumers_id = 1:100))

  dimnames(cost) <- list(names(co), names(pr))

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

test_that("non-missing 0s are handled", {
  cost <- cost_adj <- 1 / lma_simets(data.frame(x = 1:5, y = 1:5), metric = "euc", symmetrical = TRUE) - 1
  diag(cost_adj) <- 1e-6
  expect_identical(catchment_ratio(1:5, 1:5, cost = cost), catchment_ratio(1:5, 1:5, cost = cost_adj))
  nas <- sample.int(length(cost), 5)
  zeros <- sample.int(length(cost), 2)
  cost[nas] <- cost_adj[nas] <- NA
  cost[zeros] <- 0
  cost_adj[zeros] <- 1e-6
  expect_identical(catchment_ratio(1:5, 1:5, cost = cost), catchment_ratio(1:5, 1:5, cost = cost_adj))
})

test_that("verbose works", {
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

  expect_identical(
    sub("^(?:[^\\s]{1}|[^a-z]+)\\s", "", capture.output(
      catchment_ratio(supply, demand, verbose = TRUE),
      type = "message"
    ))[2:8],
    c(
      "consumers value: 1",
      "providers value: 1",
      "consumers id: sequence along consumers value",
      "providers id: sequence along providers value",
      "cost: 1",
      "weight: cost",
      "calculated 2-step floating catchment area (resources per consumer)"
    )
  )
  expect_identical(
    sub("^(?:[^\\s]{1}|[^a-z]+)\\s", "", capture.output(
      catchment_ratio(
        demand, supply,
        verbose = TRUE, normalize_weight = TRUE,
        consumers_id = "id", consumers_value = "d", consumers_location = c("lat", "long"),
        providers_id = "id", providers_value = "s", providers_location = c("lat", "long"),
      ),
      type = "message"
    ))[2:12],
    c(
      "consumers value: `d` column",
      "providers value: `s` column",
      "consumers id: `id` column",
      "providers id: `id` column",
      "calculating cost from locations...",
      "consumers location: `consumers` columns (lat and long)",
      "providers location: `providers` columns (lat and long)",
      "cost: calculated Euclidean distances",
      "weight: cost",
      "normalizing weight",
      "calculated 3-step floating catchment area (resources per consumer)"
    )
  )
})
