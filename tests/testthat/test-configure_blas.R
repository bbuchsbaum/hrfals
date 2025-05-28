context("configure_blas utility")

# configure_blas should invisibly return BLAS library name

lib_info <- withVisible(configure_blas(1))
expected <- extSoftVersion()["BLAS"]

test_that("configure_blas returns BLAS info invisibly", {
  expect_true(is.character(lib_info$value))
  expect_equal(lib_info$value, expected)
  expect_false(lib_info$visible)
})
