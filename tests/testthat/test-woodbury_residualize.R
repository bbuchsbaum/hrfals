context("woodbury residualization")

set.seed(123)

n <- 30
m <- 4
Tt <- 10
A <- matrix(rnorm(n * m), n, m)
C <- matrix(rnorm(n * Tt), n, Tt)

explicit_residualize <- function(C, A, lambda = 0) {
  AtA <- crossprod(A)
  m <- ncol(A)
  if (lambda != 0) AtA <- AtA + lambda * diag(m)
  P <- cholSolve(AtA, crossprod(A, C))
  C - A %*% P
}

V1 <- woodbury_residualize(C, A, lambda_ridge = 0.01)
V2 <- explicit_residualize(C, A, lambda = 0.01)

test_that("woodbury_residualize matches explicit projection", {
  expect_equal(max(abs(V1 - V2)), 0, tolerance = 1e-12)
  expect_equal(dim(V1), dim(C))
})
