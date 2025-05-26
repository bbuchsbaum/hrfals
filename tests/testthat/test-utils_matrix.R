library(testthat)
context("cholSolve utility")

set.seed(1)

M_small <- diag(c(1, 1e-12, 2))
b <- rnorm(3)

# when chol succeeds but diag element is below eps
res <- cholSolve(M_small, b, eps = 1e-5)
expect_equal(res, solve(M_small + 1e-5 * diag(3), b))

# when chol fails (rank deficient)
M_rankdef <- diag(c(1, 0, 2))
res2 <- cholSolve(M_rankdef, b, eps = 1e-5)
expect_equal(res2, solve(M_rankdef + 1e-5 * diag(3), b))
