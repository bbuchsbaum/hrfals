context("cf_als_engine")

simple_cfals_data <- function() {
  set.seed(123)
  n <- 50
  d <- 3
  k <- 2
  v <- 4
  h_true <- matrix(rnorm(d * v), d, v)
  beta_true <- matrix(rnorm(k * v), k, v)
  X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
  Xbig <- do.call(cbind, X_list)
  Y <- Xbig %*% as.vector(matrix(h_true, d, v) %*% t(beta_true))
  Y <- matrix(Y, n, v)
  list(X_list = X_list, Y = Y, d = d, k = k)
}

test_that("cf_als_engine returns matrices with correct dimensions", {
  dat <- simple_cfals_data()
  res <- cf_als_engine(dat$X_list, dat$Y,
                       lambda_b = 0.1,
                       lambda_h = 0.1,
                       R_mat_eff = NULL,
                       fullXtX_flag = FALSE,
                       precompute_xty_flag = TRUE,
                       max_alt = 1)
  expect_equal(dim(res$h), c(dat$d, ncol(dat$Y)))
  expect_equal(dim(res$beta), c(dat$k, ncol(dat$Y)))
})

simple_small_data <- function() {
  set.seed(42)
  n <- 20; d <- 2; k <- 2; v <- 2
  h_true <- matrix(rnorm(d * v), d, v)
  beta_true <- matrix(rnorm(k * v), k, v)
  X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
  Y <- matrix(0, n, v)
  for (c in seq_len(k)) {
    Y <- Y + (X_list[[c]] %*% h_true) *
      matrix(rep(beta_true[c, ], each = n), n, v)
  }
  list(X_list = X_list, Y = Y)
}

test_that("XtY strategies give identical results", {
  dat <- simple_small_data()
  Rm <- diag(2) * 1.5
  res_pre <- cf_als_engine(dat$X_list, dat$Y,
                           lambda_b = 0.1,
                           lambda_h = 0.2,
                           R_mat_eff = Rm,
                           fullXtX_flag = FALSE,
                           precompute_xty_flag = TRUE,
                           max_alt = 1)
  res_onfly <- cf_als_engine(dat$X_list, dat$Y,
                             lambda_b = 0.1,
                             lambda_h = 0.2,
                             R_mat_eff = Rm,
                             fullXtX_flag = FALSE,
                             precompute_xty_flag = FALSE,
                             max_alt = 1)
  expect_equal(res_pre$h, res_onfly$h, tolerance = 1e-12)
  expect_equal(res_pre$beta, res_onfly$beta, tolerance = 1e-12)
})

test_that("XtY strategies match with fullXtX", {
  dat <- simple_small_data()
  Rm <- diag(2) * 1.5
  res_pre <- cf_als_engine(dat$X_list, dat$Y,
                           lambda_b = 0.1,
                           lambda_h = 0.2,
                           R_mat_eff = Rm,
                           fullXtX_flag = TRUE,
                           precompute_xty_flag = TRUE,
                           max_alt = 1)
  res_onfly <- cf_als_engine(dat$X_list, dat$Y,
                             lambda_b = 0.1,
                             lambda_h = 0.2,
                             R_mat_eff = Rm,
                             fullXtX_flag = TRUE,
                             precompute_xty_flag = FALSE,
                             max_alt = 1)
  expect_equal(res_pre$h, res_onfly$h, tolerance = 1e-12)
  expect_equal(res_pre$beta, res_onfly$beta, tolerance = 1e-12)
})
