context("sparse CFALS stress test")

# helper function to simulate large data set
simulate_large_predictor_data <- function(n = 500, d = 20, k = 500,
                                          n_active = 10, noise_sd = 0.05) {
  set.seed(123)
  h_true <- rnorm(d)
  beta_true <- c(rep(1, n_active), rep(0, k - n_active))
  X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
  Y <- matrix(0, n, 1)
  for (c in seq_len(k)) {
    Y <- Y + (X_list[[c]] %*% h_true) * beta_true[c]
  }
  Y <- Y + matrix(rnorm(n, sd = noise_sd), n, 1)
  list(X_list = X_list, Y = Y, d = d, k = k)
}


test_that("cf_als_engine handles many predictors with caching options", {
  skip_on_cran()
  dat <- simulate_large_predictor_data()

  res_cache <- cf_als_engine(
    dat$X_list, dat$Y,
    lambda_b = 0.1,
    lambda_h = 0.1,
    lambda_init = 0.1,
    R_mat_eff = diag(dat$d),
    fullXtX_flag = FALSE,
    precompute_xty_flag = TRUE,
    Phi_recon_matrix = diag(dat$d),
    h_ref_shape_canonical = rep(1, dat$d),
    max_alt = 1,
    beta_penalty = list(l1 = 0.05, alpha = 1, warm_start = TRUE),
    design_control = list(standardize_predictors = FALSE,
                          cache_design_blocks = TRUE)
  )

  res_nocache <- cf_als_engine(
    dat$X_list, dat$Y,
    lambda_b = 0.1,
    lambda_h = 0.1,
    lambda_init = 0.1,
    R_mat_eff = diag(dat$d),
    fullXtX_flag = FALSE,
    precompute_xty_flag = TRUE,
    Phi_recon_matrix = diag(dat$d),
    h_ref_shape_canonical = rep(1, dat$d),
    max_alt = 1,
    beta_penalty = list(l1 = 0.05, alpha = 1, warm_start = TRUE),
    design_control = list(standardize_predictors = FALSE,
                          cache_design_blocks = FALSE)
  )

  expect_equal(dim(res_cache$beta), c(dat$k, 1))
  expect_equal(dim(res_nocache$beta), c(dat$k, 1))
  expect_equal(dim(res_cache$h), c(dat$d, 1))
  expect_equal(dim(res_nocache$h), c(dat$d, 1))
})
