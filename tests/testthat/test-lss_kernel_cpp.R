context("lss_kernel_cpp")

naive_lss_kernel <- function(C, A, Y, lambda_ridge = 0) {
  AtA <- crossprod(A)
  if (lambda_ridge != 0) AtA <- AtA + lambda_ridge * diag(ncol(A))
  P <- cholSolve(AtA, t(A))
  Cres <- C - A %*% (P %*% C)
  Yres <- Y - A %*% (P %*% Y)
  B <- matrix(0, ncol(C), ncol(Y))
  for (t in seq_len(ncol(C))) {
    denom <- sum(Cres[,t]^2)
    if (denom != 0) {
      num <- crossprod(Cres[,t], Yres)
      B[t,] <- num / denom
    }
  }
  B
}

simple_kernel_data <- function() {
  set.seed(1)
  n <- 20; m <- 4; Tt <- 5; v <- 3
  C <- matrix(rnorm(n*Tt), n, Tt)
  A <- matrix(rnorm(n*m), n, m)
  Y <- matrix(rnorm(n*v), n, v)
  list(C=C,A=A,Y=Y)
}

test_that("lss_kernel_cpp matches naive implementation", {
  dat <- simple_kernel_data()
  res_cpp <- lss_kernel_cpp(dat$C, dat$A, dat$Y, lambda_ridge = 0.1)
  res_naive <- naive_lss_kernel(dat$C, dat$A, dat$Y, lambda_ridge = 0.1)
  expect_equal(res_cpp, res_naive, tolerance = 1e-12)
})
