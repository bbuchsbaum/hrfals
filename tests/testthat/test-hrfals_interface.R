context("hrfals interface")

library(fmridesign)

test_that("hrfals_from_design produces valid hrfals_fit object", {
  dat <- simulate_cfals_wrapper_data(fmrihrf::HRF_SPMG3)
  design <- create_cfals_design(dat$Y, dat$event_model, fmrihrf::HRF_SPMG3)

  ctrl <- list(lambda_init = 1, lambda_b = 0.1, lambda_h = 0.1,
               fullXtX = TRUE, max_alt = 2)

  # Test the hrfals_from_design interface function
  fit <- hrfals_from_design(dat$Y, design, method = "cf_als", control = ctrl)
  
  # Check that the fit object has all expected components
  expect_s3_class(fit, "hrfals_fit")
  expect_equal(dim(fit$h_coeffs), c(fmrihrf::nbasis(fmrihrf::HRF_SPMG3), ncol(dat$Y)))
  expect_equal(dim(fit$beta_amps), c(length(design$X_list), ncol(dat$Y)))
  expect_equal(fit$lambdas[["beta"]], 0.1)
  expect_equal(fit$lambdas[["h"]], 0.1)
  expect_equal(fit$design_info$fullXtX, TRUE)
  
  # Check that the fit produces reasonable results
  expect_true(all(is.finite(fit$h_coeffs)))
  expect_true(all(is.finite(fit$beta_amps)))
  gof_finite <- fit$gof[is.finite(fit$gof)]
  expect_true(all(gof_finite <= 1))  # R² cannot exceed 1, but can be arbitrarily negative
})

test_that("hrfals_from_design forwards lambda_init to cf_als_engine", {
  dat <- simple_small_data()
  design <- list(
    X_list_proj = dat$X_list,
    Y_proj = dat$Y,
    Phi_recon_matrix = dat$Phi,
    h_ref_shape_canonical = dat$href,
    d_basis_dim = ncol(dat$X_list[[1]]),
    k_conditions = length(dat$X_list),
    n_timepoints = nrow(dat$Y),
    v_voxels = ncol(dat$Y),
    predictor_means = rep(0, length(dat$X_list)),
    predictor_sds = rep(1, length(dat$X_list)),
    hrf_basis = fmrihrf::HRF_SPMG2
  )

  fit_lo <- hrfals_from_design(dat$Y, design, method = "cf_als",
                               control = list(lambda_init = 0.001, max_alt = 1,
                                              lambda_b = 0.1, lambda_h = 0.1))
  fit_hi <- hrfals_from_design(dat$Y, design, method = "cf_als",
                               control = list(lambda_init = 100, max_alt = 1,
                                              lambda_b = 0.1, lambda_h = 0.1))

  expect_gt(max(abs(fit_lo$h_coeffs - fit_hi$h_coeffs)), 1e-10)
})

test_that("hrfals_from_design rescales betas using design$predictor_sds", {
  set.seed(123)
  n <- 30; d <- 2; k <- 2; v <- 2
  scales <- c(cond1 = 2, cond2 = 3)
  X_orig <- list(
    cond1 = matrix(rnorm(n * d), n, d),
    cond2 = matrix(rnorm(n * d), n, d)
  )
  X_std <- list(
    cond1 = X_orig$cond1 / scales[["cond1"]],
    cond2 = X_orig$cond2 / scales[["cond2"]]
  )
  h_true <- matrix(rnorm(d * v), d, v)
  beta_orig <- matrix(rnorm(k * v), k, v)
  beta_std <- sweep(beta_orig, 1, scales, "*")

  Y <- matrix(0, n, v)
  for (c in seq_len(k)) {
    Y <- Y + (X_std[[c]] %*% h_true) *
      matrix(rep(beta_std[c, ], each = n), n, v)
  }

  Phi <- diag(d)
  href <- rep(1, d)
  res_engine <- cf_als_engine(
    X_std, Y,
    lambda_init = 1,
    lambda_b = 0.1,
    lambda_h = 0.1,
    Phi_recon_matrix = Phi,
    h_ref_shape_canonical = href,
    max_alt = 1,
    design_control = list(standardize_predictors = FALSE, cache_design_blocks = TRUE)
  )

  design <- list(
    X_list_proj = X_std,
    Y_proj = Y,
    Phi_recon_matrix = Phi,
    h_ref_shape_canonical = href,
    d_basis_dim = d,
    k_conditions = k,
    n_timepoints = n,
    v_voxels = v,
    predictor_means = c(cond1 = 0, cond2 = 0),
    predictor_sds = scales,
    hrf_basis = fmrihrf::HRF_SPMG2,
    condition_names = names(X_std)
  )

  fit <- hrfals_from_design(Y, design, method = "cf_als",
                            control = list(lambda_init = 1, max_alt = 1,
                                           lambda_b = 0.1, lambda_h = 0.1))

  expected_beta <- sweep(res_engine$beta, 1, scales, "/")
  expect_equal(fit$beta_amps, expected_beta)
})
