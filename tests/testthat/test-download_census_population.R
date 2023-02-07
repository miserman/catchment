test_that("files are downloaded and loaded", {
  dir <- "../../census_download_tests"
  if (grepl("R_LIBS", getwd(), fixed = TRUE)) dir.create(dir, FALSE, TRUE)
  if (!dir.exists(dir)) dir <- "../../../../census_download_tests"
  skip_if_not(dir.exists(dir), "not downloading data")
  files <- paste0(dir, "/", c("dc_population", "dc_population_margins", "dc_commutes"), "_2019.csv")
  unlink(files)
  data <- download_census_population(dir, "dc", include_margins = TRUE, include_commutes = TRUE)
  names(files) <- c("estimates", "margins", "commutes")
  expect_identical(names(data), names(files))
  expect_true(all(file.exists(files)))
  for (f in names(files)) {
    rd <- read.csv(files[[f]], row.names = 1)
    if (colnames(data[[f]])[1] == "GEOID") {
      rownames(data[[f]]) <- data[[f]][, 1]
      data[[f]] <- data[[f]][, -1]
    }
    expect_identical(rownames(rd), rownames(data[[f]]))
    expect_identical(as.numeric(as.matrix(rd)), as.numeric(as.matrix(data[[f]])))
  }
})
