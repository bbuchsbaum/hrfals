context("estimate_hrf_cfals wrapper")

library(fmrireg)

simulate_cfals_wrapper_data <- function(hrf_basis, noise_sd = 0.05, signal_scale = 1) {
  sf <- sampling_frame(blocklens = 60, TR = 1)
  events <- data.frame(
    onset = c(5, 15, 30, 45),
    condition = factor(c("A", "A", "B", "B")),
    block = 1
  )
  emod <- event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  
  # Use create_fmri_design to properly create design matrices
  design <- create_fmri_design(emod, hrf_basis)
  X_list <- design$X_list
  
  d <- design$d
  k <- design$k
  v <- 2
  n_timepoints <- length(samples(sf, global = TRUE))
  
  h_true <- matrix(rnorm(d * v), d, v) * signal_scale
  beta_true <- matrix(rnorm(k * v), k, v) * signal_scale
  Y <- matrix(0, n_timepoints, v)
  for (c in seq_along(X_list)) {
    Y <- Y + (X_list[[c]] %*% h_true) *
      matrix(rep(beta_true[c, ], each = nrow(Y)), nrow(Y), v)
  }
  Y <- Y + matrix(rnorm(length(Y), sd = noise_sd), n_timepoints, v)
  attr(Y, "sampling_frame") <- sf
  list(Y = Y, event_model = emod, X_list = X_list,
       h_true = h_true, beta_true = beta_true, sframe = sf)
}


test_that("estimate_hrf_cfals returns expected dimensions", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  fit <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", HRF_SPMG3,
                            lambda_b = 0.1, lambda_h = 0.1)
  expect_s3_class(fit, "hrfals_fit")
  expect_equal(dim(fit$h_coeffs), c(nbasis(HRF_SPMG3), ncol(dat$Y)))
  expect_equal(dim(fit$beta_amps), c(2, ncol(dat$Y)))
  expect_equal(rownames(fit$beta_amps), c("conditionA", "conditionB"))
  expect_equal(fit$target_event_term_name, "hrf(condition)")
  expect_true(is.matrix(fit$phi_recon_matrix))
})

test_that("estimate_hrf_cfals carries bad_row_idx", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  dat$Y[4, 1] <- NA
  fit <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", HRF_SPMG3,
                            lambda_b = 0.1, lambda_h = 0.1)
  expect_equal(fit$bad_row_idx, 4)
})


test_that("estimate_hrf_cfals matches direct ls_svd_1als", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  prep <- create_cfals_design(dat$Y, dat$event_model, HRF_SPMG3)
  direct <- ls_svd_1als_engine(prep$X_list_proj, prep$Y_proj,
                               lambda_init = 0,
                               lambda_b = 0.1,
                               lambda_h = 0.1,
                               fullXtX_flag = TRUE,
                               Phi_recon_matrix = prep$Phi_recon_matrix,
                               h_ref_shape_canonical = prep$h_ref_shape_canonical,
                               R_mat = diag(prep$d_basis_dim))
  wrap <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", HRF_SPMG3,
                             method = "ls_svd_1als",
                             lambda_init = 0,
                             lambda_b = 0.1,
                             lambda_h = 0.1,
                             fullXtX = TRUE,
                             R_mat = "identity")
  expect_equal(wrap$h_coeffs, direct$h)
  expect_equal(wrap$beta_amps, direct$beta)
})



test_that("estimate_hrf_cfals predictions match canonical GLM", {
  set.seed(123)
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  fit <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)",
                            HRF_SPMG3,
                            method = "cf_als",
                            lambda_b = 0,
                            lambda_h = 0,
                            max_alt = 1)
  n <- nrow(dat$Y)
  v <- ncol(dat$Y)
  pred_cfals <- matrix(0, n, v)
  for (c in seq_along(dat$X_list)) {
    pred_cfals <- pred_cfals + (dat$X_list[[c]] %*% fit$h_coeffs) *
      matrix(rep(fit$beta_amps[c, ], each = n), n, v)
  }
  Xbig <- do.call(cbind, dat$X_list)
  gamma_hat <- chol2inv(chol(crossprod(Xbig))) %*% crossprod(Xbig, dat$Y)
  pred_glm <- Xbig %*% gamma_hat
  expect_equal(pred_cfals, pred_glm, tolerance = 5e-2)
})
                   
test_that("R_mat = 'basis_default' uses basis penalty matrix", {

  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  prep <- create_cfals_design(dat$Y, dat$event_model, HRF_SPMG3)
  Rb <- fmrireg::penalty_matrix(HRF_SPMG3)
  direct <- ls_svd_1als_engine(prep$X_list_proj, prep$Y_proj,
                               lambda_init = 0,
                               lambda_b = 0.1,
                               lambda_h = 0.1,
                               fullXtX_flag = TRUE,
                               Phi_recon_matrix = prep$Phi_recon_matrix,
                               h_ref_shape_canonical = prep$h_ref_shape_canonical,
                               R_mat = Rb)
  wrap <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", HRF_SPMG3,
                             method = "ls_svd_1als",
                             lambda_init = 0,
                             lambda_b = 0.1,
                             lambda_h = 0.1,
                             fullXtX = TRUE,
                             R_mat = "basis_default")
  expect_equal(wrap$h_coeffs, direct$h)
  expect_equal(wrap$beta_amps, direct$beta)
})

test_that("R_mat custom matrix is used", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  prep <- create_cfals_design(dat$Y, dat$event_model, HRF_SPMG3)
  R_custom <- diag(prep$d_basis_dim) * 2
  direct <- ls_svd_1als_engine(prep$X_list_proj, prep$Y_proj,
                               lambda_init = 0,
                               lambda_b = 0.1,
                               lambda_h = 0.1,
                               fullXtX_flag = TRUE,
                               Phi_recon_matrix = prep$Phi_recon_matrix,
                               h_ref_shape_canonical = prep$h_ref_shape_canonical,
                               R_mat = R_custom)
  wrap <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", HRF_SPMG3,
                             method = "ls_svd_1als",
                             lambda_init = 0,
                             lambda_b = 0.1,
                             lambda_h = 0.1,
                             fullXtX = TRUE,
                             R_mat = R_custom)
  expect_equal(wrap$h_coeffs, direct$h)
  expect_equal(wrap$beta_amps, direct$beta)

})


simulate_multiterm_data <- function(hrf_basis, noise_sd = 0.05) {
  sf <- sampling_frame(blocklens = 60, TR = 1)
  events <- data.frame(
    onset = c(5, 15, 25, 35),
    term1 = factor(c("A", "A", "B", "B")),
    term2 = factor(c("C", "D", "C", "D")),
    block = 1
  )
  emod <- event_model(onset ~ hrf(term1) + hrf(term2), data = events,
                      block = ~ block, sampling_frame = sf)
  
  # Use create_fmri_design to properly create design matrices
  design <- create_fmri_design(emod, hrf_basis)
  X_list <- design$X_list
  
  d <- design$d
  k <- design$k
  v <- 2
  n_timepoints <- length(samples(sf, global = TRUE))
  
  h_true <- matrix(rnorm(d * v), d, v)
  beta_true <- matrix(rnorm(k * v), k, v)
  Y <- matrix(0, n_timepoints, v)
  for (c in seq_along(X_list)) {
    Y <- Y + (X_list[[c]] %*% h_true) *
      matrix(rep(beta_true[c, ], each = nrow(Y)), nrow(Y), v)
  }
  Y <- Y + matrix(rnorm(length(Y), sd = noise_sd), n_timepoints, v)
  attr(Y, "sampling_frame") <- sf
  list(Y = Y, event_model = emod)
}


test_that("estimate_hrf_cfals integrates across HRF bases and terms", {
  bases <- list(HRF_SPMG3, gen_hrf(hrf_bspline, N=4))
  for (b in bases) {
    dat <- simulate_multiterm_data(b)
    for (term in c("hrf(term1)", "hrf(term2)")) {
      fit <- estimate_hrf_cfals(dat$Y, dat$event_model, term, b,
                                lambda_b = 0.1, lambda_h = 0.1)
      expect_s3_class(fit, "hrfals_fit")
      expect_equal(nrow(fit$h_coeffs), nbasis(b))
      expect_equal(fit$target_event_term_name, term)
    }
  }
})

test_that("R_mat options work", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)

  fit_basis <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)",
                                  HRF_SPMG3,
                                  lambda_b = 0.1, lambda_h = 0.1,
                                  R_mat = "basis_default")
  expect_s3_class(fit_basis, "hrfals_fit")

  Rm <- diag(nbasis(HRF_SPMG3))
  fit_custom <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)",
                                   HRF_SPMG3,
                                   lambda_b = 0.1, lambda_h = 0.1,
                                   R_mat = Rm)
  expect_s3_class(fit_custom, "hrfals_fit")

  expect_error(
    estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", HRF_SPMG3,
                       lambda_b = 0.1, lambda_h = 0.1,
                       R_mat = NA),
    "R_mat must be 'identity', 'basis_default', or a numeric matrix"
  )
})

test_that("estimate_hrf_spatial_cfals forwards arguments", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  mask <- array(1, dim = c(2, 1, 1))
  lap_obj <- build_voxel_laplacian(mask)
  fit <- estimate_hrf_spatial_cfals(dat$Y, dat$event_model,
                                    "hrf(condition)", HRF_SPMG3,
                                    laplacian_obj = lap_obj,
                                    lambda_s_default = 0.05,
                                    h_solver = "direct",
                                    max_alt = 1)
  expect_s3_class(fit, "hrfals_fit")
})

