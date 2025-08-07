load_auto_read_csv <- function() {
  path <- testthat::test_path("..", "..", "make_reports.R")
  exprs <- parse(path)
  env <- new.env(parent = globalenv())
  env$str_count <- stringr::str_count
  for (expr in exprs) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("auto_read_csv"))) {
      eval(expr, envir = env)
      return(env$auto_read_csv)
    }
  }
  stop("auto_read_csv not found")
}

auto_read_csv <- load_auto_read_csv()

testthat::test_that("reads comma-separated CSV", {
  path <- testthat::test_path("fixtures", "comma.csv")
  df <- auto_read_csv(path)
  testthat::expect_equal(colnames(df), c("a", "b", "c"))
})

testthat::test_that("reads semicolon-separated CSV", {
  path <- testthat::test_path("fixtures", "semicolon.csv")
  df <- auto_read_csv(path)
  testthat::expect_equal(colnames(df), c("a", "b", "c"))
})
