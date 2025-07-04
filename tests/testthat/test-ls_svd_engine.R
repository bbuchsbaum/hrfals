library(testthat)

context("ls_svd_engine")

simple_ls_svd_data <- function() {
  set.seed(123)
  n <- 50
  d <- 3
  k <- 2
  v <- 5
  h_true <- matrix(rnorm(d * v), d, v)
  beta_true <- matrix(rnorm(k * v), k, v)
  X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
  Xbig <- do.call(cbind, X_list)
  Y <- Xbig %*% as.vector(matrix(h_true, d, v) %*% t(beta_true))
  Y <- matrix(Y, n, v)
  phi <- diag(d)
  href <- rep(1, nrow(phi))
  list(X_list = X_list, Y = Y, d = d, k = k,
       Phi = phi, href = href)
}


test_that("ls_svd_engine returns matrices with correct dimensions", {
  dat <- simple_ls_svd_data()
  res <- ls_svd_engine(dat$X_list, dat$Y,
                       lambda_init = 0,
                       Phi_recon_matrix = dat$Phi,
                       h_ref_shape_canonical = dat$href)
  expect_equal(dim(res$h), c(dat$d, ncol(dat$Y)))
  expect_equal(dim(res$beta), c(dat$k, ncol(dat$Y)))
  expect_equal(dim(res$Gamma_hat), c(dat$d * dat$k, ncol(dat$Y)))
})

test_that("ls_svd_engine applies simple ridge", {
  dat <- simple_ls_svd_data()
  res <- ls_svd_engine(dat$X_list, dat$Y,
                       lambda_init = 0.5,
                       Phi_recon_matrix = dat$Phi,
                       h_ref_shape_canonical = dat$href)
  Xbig <- do.call(cbind, dat$X_list)
  XtX <- crossprod(Xbig)
  Xty <- crossprod(Xbig, dat$Y)
  Gamma_manual <- solve(XtX + 0.5 * diag(dat$d * dat$k), Xty)
  expect_equal(res$Gamma_hat, Gamma_manual)
})

test_that("ls_svd_engine requires normalised h_ref", {
  dat <- simple_ls_svd_data()
  bad_ref <- dat$href * 2
  expect_error(
    ls_svd_engine(dat$X_list, dat$Y,
                  lambda_init = 0,
                  Phi_recon_matrix = dat$Phi,
                  h_ref_shape_canonical = bad_ref),
    "must be normalised"
  )
})
