context("qr_residualize function")

set.seed(456)

n <- 30
m <- 6
Tt <- 8
A <- matrix(rnorm(n * m), n, m)
C <- matrix(rnorm(n * Tt), n, Tt)

# reference residualization using base qr.resid
qr_obj <- qr(A, LAPACK = FALSE)
ref <- qr.resid(qr_obj, C)

res_default <- qr_residualize(C, A)
res_lapack <- qr_residualize(C, A, lapack_qr = TRUE)

test_that("qr_residualize matches qr.resid", {
  expect_equal(res_default, ref, tolerance = 1e-12)
  expect_equal(dim(res_default), dim(C))
})

test_that("lapack_qr parameter gives identical result", {
  expect_equal(res_lapack, ref, tolerance = 1e-12)
})
