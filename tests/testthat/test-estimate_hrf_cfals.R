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
  reg_lists <- lapply(emod$terms, regressors.event_term,
                      hrf = hrf_basis,
                      sampling_frame = sf,
                      summate = FALSE,
                      drop.empty = TRUE)
  regs <- unlist(reg_lists, recursive = FALSE)
  sample_times <- samples(sf, global = TRUE)
  X_list <- lapply(regs, function(r)
    evaluate(r, sample_times, precision = sf$precision))
  d <- nbasis(hrf_basis)
  k <- length(X_list)
  v <- 2
  h_true <- matrix(rnorm(d * v), d, v) * signal_scale
  beta_true <- matrix(rnorm(k * v), k, v) * signal_scale
  Y <- matrix(0, nrow(sample_times), v)
  for (c in seq_along(X_list)) {
    Y <- Y + (X_list[[c]] %*% h_true) *
      matrix(rep(beta_true[c, ], each = nrow(Y)), nrow(Y), v)
  }
  Y <- Y + matrix(rnorm(length(Y), sd = noise_sd), nrow(Y), v)
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
##<<<<<<< codex/update-design-object-and-engine-arguments
                               h_ref_shape_norm = NULL,
##<<<<<<< codex/update-unit-and-wrapper-tests
##=======
                               R_mat = diag(prep$d_basis_dim),
                               Phi_recon_matrix = prep$Phi_recon_matrix,
                               h_ref_shape_canonical = prep$h_ref_shape_canonical)
##=======
##>>>>>>> main
                               Phi_recon_matrix = prep$Phi_recon_matrix,
                               h_ref_shape_canonical = prep$h_ref_shape_canonical,
                               R_mat = diag(prep$d_basis_dim))
##>>>>>>> main
  wrap <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", HRF_SPMG3,
                             method = "ls_svd_1als",
                             lambda_init = 0,
                             lambda_b = 0.1,
                             lambda_h = 0.1,
                             fullXtX = TRUE,
                             penalty_R_mat_type = "identity")
  expect_equal(wrap$h_coeffs, direct$h)
  expect_equal(wrap$beta_amps, direct$beta)
})

test_that("penalty_R_mat_type 'basis_default' uses basis penalty matrix", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  prep <- create_cfals_design(dat$Y, dat$event_model, HRF_SPMG3)
  Rb <- penalty_matrix(HRF_SPMG3)
  direct <- ls_svd_1als_engine(prep$X_list_proj, prep$Y_proj,
                               lambda_init = 0,
                               lambda_b = 0.1,
                               lambda_h = 0.1,
                               fullXtX_flag = TRUE,
                               h_ref_shape_norm = NULL,
                               R_mat = Rb)
  wrap <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", HRF_SPMG3,
                             method = "ls_svd_1als",
                             lambda_init = 0,
                             lambda_b = 0.1,
                             lambda_h = 0.1,
                             fullXtX = TRUE,
                             penalty_R_mat_type = "basis_default")
  expect_equal(wrap$h_coeffs, direct$h)
  expect_equal(wrap$beta_amps, direct$beta)
})

test_that("penalty_R_mat_type 'custom' uses provided matrix", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  prep <- create_cfals_design(dat$Y, dat$event_model, HRF_SPMG3)
  R_custom <- diag(prep$d_basis_dim) * 2
  direct <- ls_svd_1als_engine(prep$X_list_proj, prep$Y_proj,
                               lambda_init = 0,
                               lambda_b = 0.1,
                               lambda_h = 0.1,
                               fullXtX_flag = TRUE,
                               h_ref_shape_norm = NULL,
                               R_mat = R_custom)
  wrap <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", HRF_SPMG3,
                             method = "ls_svd_1als",
                             lambda_init = 0,
                             lambda_b = 0.1,
                             lambda_h = 0.1,
                             fullXtX = TRUE,
                             penalty_R_mat_type = "custom",
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
  reg_lists <- lapply(emod$terms, regressors.event_term,
                      hrf = hrf_basis,
                      sampling_frame = sf,
                      summate = FALSE,
                      drop.empty = TRUE)
  regs <- unlist(reg_lists, recursive = FALSE)
  sample_times <- samples(sf, global = TRUE)
  X_list <- lapply(regs, function(r)
    evaluate(r, sample_times, precision = sf$precision))
  d <- nbasis(hrf_basis)
  k <- length(X_list)
  v <- 2
  h_true <- matrix(rnorm(d * v), d, v)
  beta_true <- matrix(rnorm(k * v), k, v)
  Y <- matrix(0, nrow(sample_times), v)
  for (c in seq_along(X_list)) {
    Y <- Y + (X_list[[c]] %*% h_true) *
      matrix(rep(beta_true[c, ], each = nrow(Y)), nrow(Y), v)
  }
  Y <- Y + matrix(rnorm(length(Y), sd = noise_sd), nrow(Y), v)
  attr(Y, "sampling_frame") <- sf
  list(Y = Y, event_model = emod)
}


test_that("estimate_hrf_cfals integrates across HRF bases and terms", {
  bases <- list(HRF_SPMG3, hrfspline_generator(nbasis = 4))
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

