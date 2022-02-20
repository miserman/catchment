test_that("files are downloaded and loaded", {
  dir <- "../../census_download_tests"
  if (grepl("R_LIBS", getwd(), fixed = TRUE)) dir.create(dir, FALSE, TRUE)
  if (!dir.exists(dir)) dir <- "../../../../census_download_tests"
  skip_if_not(dir.exists(dir), "not downloading data")
  library(sf)
  data <- download_census_shapes(
    dir, "dc", "county",
    strip_features = TRUE,
    simplify = st_simplify, preserveTopology = TRUE
  )
  expect_identical(names(data), c("GEOID", "geometry"))
  file <- paste0(dir, "/cb_2019_us_county_500k.geojson")
  expect_true(file.exists(file))
  expect_identical(data$GEOID, st_read(file, quiet = TRUE)$GEOID)
})
