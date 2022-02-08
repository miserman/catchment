test_that("weights align with manual", {
  cost <- matrix(abs(rnorm(100 * 20, 2)), 100)
  expect_true(all((cost > 0 & cost < 1) * 1 == catchment_weight(cost, 1)))
  w <- (cost > 0 & cost <= 1) * .5
  w[cost > 0 & cost <= .5] <- 1
  expect_true(all(w == catchment_weight(cost, list(c(1, .5), c(.5, 1)))))
  expect_true(all(exp(-cost * 2) == catchment_weight(cost, "exponential")))
  expect_true(all(dnorm(cost, 0, 2) == catchment_weight(cost, "dnorm")))
  expect_true(all(dnorm(cost) == catchment_weight(cost, dnorm)))
  expect_true(all(pnorm(-cost, 0, 2) == catchment_weight(cost, "pnorm")))
  expect_true(all(pnorm(-cost) == catchment_weight(cost, function(x) pnorm(-x))))
  expect_true(all(sqrt(1 / cost^2) == catchment_weight(cost, "gravity")))
  expect_true(all(sqrt(1 / cost^2) == catchment_weight(cost, "normal")))
  expect_true(all(1 / (1 + log2(cost)) == catchment_weight(cost, "logarithmic")))
  expect_true(all({
    tw <- (2 - cost) / 2
    tw[tw < 0] <- 0
    tw
  } == catchment_weight(cost, "linear")))
  expect_true(all(1 / (1 + exp(cost * 2)) == catchment_weight(cost, "logistic")))
  expect_true(all(exp(-cost^2 / 8) == catchment_weight(cost, "gaussian")))

  w <- (cost > 0 & cost <= 1) * .5
  w[cost > 0 & cost <= .5] <- 1
  wr <- rowSums(w)
  wr[wr == 0] <- 1
  expect_true(all(w == catchment_weight(cost, list(c(1, .5), c(.5, 1)))))

  wr <- rowSums(w)
  wr[wr == 0] <- 1
  expect_true(all(
    w * sweep(w, 1, wr, "/") ==
    catchment_weight(cost, list(c(1, .5), c(.5, 1)), normalize_weight = TRUE)
  ))

  w[cost > 3] <- 0
  expect_true(all(w == catchment_weight(cost, list(c(1, .5), c(.5, 1)), max_cost = 3)))
})

test_that("alternate inputs are handled", {
  x <- abs(rnorm(100))
  wg <- sqrt(1 / x ^ 2)
  expect_identical(wg, as.numeric(catchment_weight(x, "gravity")))
  expect_identical(
    wg,
    as.numeric(catchment_weight(as.data.frame(matrix(x, 10)), "gravity"))
  )
})

test_that("NA inputs are handled", {
  x <- abs(rnorm(100))
  x[sample.int(100, 10)] <- NA
  xz <- x
  xz[is.na(x)] <- 0
  expect_identical(catchment_weight(x, "linear"), catchment_weight(xz, "linear"))
})
