test_that("results align with manual", {
  higher <- data.frame(id = c(100, 110), pop = c(500, 500))
  lower <- data.frame(id = paste0(higher$id, 1:5), value = 1:5 / 10, pop = c(50, 10, 100, 40, 30))
  expect_identical(
    c(
      "100" = sum(lower[c(1, 3, 5), "value"] * lower[c(1, 3, 5), "pop"]) / higher[1, "pop"],
      "110" = sum(lower[c(2, 4), "value"] * lower[c(2, 4), "pop"]) / higher[2, "pop"]
    ),
    catchment_aggregate(lower, higher, consumers = "pop", id = "id", value = "value", to_id = "id")
  )
  expect_identical(
    c(
      "100" = sum(lower[c(1, 3, 5), "value"] * lower[c(1, 3, 5), "pop"]),
      "110" = sum(lower[c(2, 4), "value"] * lower[c(2, 4), "pop"])
    ),
    catchment_aggregate(
      lower, higher,
      consumers = "pop", id = "id", value = "value", to_id = "id",
      return_type = "region"
    )
  )
  expect_identical(
    c(
      "100" = sum(lower[c(1, 3, 5), "value"] * lower[c(1, 3, 5), "pop"]) / sum(lower[c(1, 3, 5), "pop"]),
      "110" = sum(lower[c(2, 4), "value"] * lower[c(2, 4), "pop"]) / sum(lower[c(2, 4), "pop"])
    ),
    catchment_aggregate(lower, 3, consumers = "pop", id = "id", value = "value", verbose = TRUE)
  )
  expect_identical(
    c("100" = sum(lower[c(1, 3, 5), "value"]), "110" = sum(lower[c(2, 4), "value"])),
    catchment_aggregate(
      lower, higher,
      consumers = "pop", id = "id", value = "value", to_id = "id",
      to_consumers = NULL, original_from = FALSE, return_type = "region"
    )
  )
  expect_identical(
    c(
      "100" = sum(lower[1:2, "value"] * lower[1:2, "pop"]) / higher[1, "pop"],
      "110" = sum(lower[3:5, "value"] * lower[3:5, "pop"]) / higher[2, "pop"]
    ),
    catchment_aggregate(
      lower, higher,
      consumers = "pop", id = "id", value = "value", to_id = "id",
      map = list("100" = lower[1:2, "id"], "110" = lower[3:5, "id"])
    )
  )
})

test_that("alternate inputs work", {
  higher <- data.frame(id = c(100, 110), pop = c(500, 500))
  lower <- data.frame(id = paste0(higher$id, 1:5), value = 1:5 / 10, pop = c(50, 10, 100, 40, 30))
  m <- c(
    "100" = sum(lower[c(1, 3, 5), "value"] * lower[c(1, 3, 5), "pop"]) / higher[1, "pop"],
    "110" = sum(lower[c(2, 4), "value"] * lower[c(2, 4), "pop"]) / higher[2, "pop"]
  )
  expect_identical(
    m * 1000,
    catchment_aggregate(
      lower, higher,
      consumers = "pop", id = "id", value = "value", to_id = "id", return_type = 1000
    )
  )
  lower$hid <- substr(lower$id, 1, 3)
  expect_identical(
    m,
    catchment_aggregate(lower, "hid", consumers = "pop", id = "id", value = "value", to_consumers = higher$pop)
  )
  expect_identical(
    m,
    catchment_aggregate(lower, higher$id, consumers = "pop", id = "id", value = "value", to_consumers = higher$pop)
  )
  v <- structure(lower$value, names = lower$id)
  expect_identical(
    m,
    catchment_aggregate(v, higher, consumers = lower$pop, to_id = "id", to_consumers = higher$pop)
  )
  expect_identical(
    m,
    catchment_aggregate(v, function(v) {
      substr(v, 1, 3)
    }, consumers = lower$pop, to_id = "id", to_consumers = higher$pop)
  )
})

test_that("verbose works", {
  higher <- data.frame(id = c(100, 110))
  lower <- data.frame(id = paste0(higher$id, 1:5), value = 1)
  expect_identical(
    sub("^(?:[^\\s]{1}|[^a-z]+)\\s", "", capture.output(
      catchment_aggregate(id = lower$id, value = lower$value, to_id = higher$id, verbose = TRUE),
      type = "message"
    )[2:6]),
    c(
      "from IDs: `id` vector",
      "from values: `value` vector",
      "to IDs: `to_id` vector",
      "mapping: by first 3 character match",
      "returning sum of value over total `from` consumers per `to` ID"
    )
  )
  expect_identical(
    sub("^(?:[^\\s]{1}|[^a-z]+)\\s", "", capture.output(
      catchment_aggregate(lower, higher, "id", verbose = TRUE, return_type = "r"),
      type = "message"
    )[2:6]),
    c(
      "from IDs: `id` column",
      "from values: 1",
      "to IDs: `id` column",
      "mapping: by first 3 character match",
      "returning sum of value per `to` ID"
    )
  )
})
