# rebuild
styler::style_pkg(filetype = c("R", "Rmd"))
spelling::spell_check_package()
devtools::document()
pkgdown::build_site()
covr::report(covr::package_coverage(quiet = FALSE), "docs/coverage.html")

# check
testthat::test_local()
devtools::check()
