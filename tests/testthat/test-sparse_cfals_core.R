context("sparse CFALS core functionality")

library(fmrireg)

# helper to simulate data with many continuous predictors and sparse betas
simulate_sparse_predictor_data <- function(n = 40, d = 3, k = 5,
                                           n_active = 2, noise_sd = 0.05) {
  set.seed(1)
  h_true <- rnorm(d)
  beta_true <- c(rep(1, n_active), rep(0, k - n_active))
  X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
  Y <- matrix(0, n, 1)
  for (c in seq_len(k)) {
    Y <- Y + (X_list[[c]] %*% h_true) * beta_true[c]
  }
  Y <- Y + matrix(rnorm(n, sd = noise_sd), n, 1)
  list(X_list = X_list, Y = Y, h_true = h_true, beta_true = beta_true)
}


test_that("l1 = 0 reproduces original results", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  fit_default <- hrfals(dat$Y, dat$event_model, HRF_SPMG3,
                        lam_beta = 0.1, lam_h = 0.1, max_alt = 1)
  fit_explicit <- hrfals(dat$Y, dat$event_model, HRF_SPMG3,
                         lam_beta = 0.1, lam_h = 0.1, max_alt = 1,
                         beta_penalty = list(l1 = 0, alpha = 1, warm_start = TRUE))
  expect_equal(fit_default$h_coeffs, fit_explicit$h_coeffs)
  expect_equal(fit_default$beta_amps, fit_explicit$beta_amps)
})


test_that("sparse beta step selects active predictors", {
  dat <- simulate_sparse_predictor_data()
  res <- cf_als_engine(dat$X_list, dat$Y,
                       lambda_b = 0.1,
                       lambda_h = 0.1,
                       lambda_init = 0.1,
                       R_mat_eff = diag(3),
                       fullXtX_flag = FALSE,
                       precompute_xty_flag = TRUE,
                       Phi_recon_matrix = diag(3),
                       h_ref_shape_canonical = rep(1, 3),
                       max_alt = 2,
                       beta_penalty = list(l1 = 0.05, alpha = 1, warm_start = TRUE),
                       design_control = list(standardize_predictors = FALSE,
                                             cache_design_blocks = TRUE))
  inactive_idx <- 3:5
  active_idx <- 1:2
  expect_true(all(abs(res$beta[inactive_idx, 1]) < 1e-6))
  expect_true(all(abs(res$beta[active_idx, 1]) > 0.1))
})


test_that("predictor standardization rescales betas", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  fit_raw <- hrfals(dat$Y, dat$event_model, HRF_SPMG3,
                    lam_beta = 0.1, lam_h = 0.1, max_alt = 1,
                    design_control = list(standardize_predictors = FALSE,
                                          cache_design_blocks = TRUE))
  fit_std <- hrfals(dat$Y, dat$event_model, HRF_SPMG3,
                    lam_beta = 0.1, lam_h = 0.1, max_alt = 1,
                    design_control = list(standardize_predictors = TRUE,
                                          cache_design_blocks = TRUE))
  expect_equal(fit_raw$beta_amps, fit_std$beta_amps, tolerance = 0.2)
})


test_that("warm_start option converges to similar solution", {
  dat <- simulate_sparse_predictor_data()
  res_ws <- cf_als_engine(dat$X_list, dat$Y,
                          lambda_b = 0.1, lambda_h = 0.1,
                          lambda_init = 0.1,
                          R_mat_eff = diag(3),
                          fullXtX_flag = FALSE,
                          precompute_xty_flag = TRUE,
                          Phi_recon_matrix = diag(3),
                          h_ref_shape_canonical = rep(1, 3),
                          max_alt = 2,
                          beta_penalty = list(l1 = 0.05, alpha = 1, warm_start = TRUE),
                          design_control = list(standardize_predictors = FALSE,
                                                cache_design_blocks = TRUE))
  res_no_ws <- cf_als_engine(dat$X_list, dat$Y,
                             lambda_b = 0.1, lambda_h = 0.1,
                             lambda_init = 0.1,
                             R_mat_eff = diag(3),
                             fullXtX_flag = FALSE,
                             precompute_xty_flag = TRUE,
                             Phi_recon_matrix = diag(3),
                             h_ref_shape_canonical = rep(1, 3),
                             max_alt = 2,
                             beta_penalty = list(l1 = 0.05, alpha = 1, warm_start = FALSE),
                             design_control = list(standardize_predictors = FALSE,
                                                   cache_design_blocks = TRUE))
  diff <- mean(abs(res_ws$beta - res_no_ws$beta))
  expect_lt(diff, 0.5)
})

