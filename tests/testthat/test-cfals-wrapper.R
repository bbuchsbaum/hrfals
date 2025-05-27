context("cfals wrapper")

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
  Y <- Y + matrix(rnorm(length(Y), sd = noise_sd), nrow(Y), v)
  attr(Y, "sampling_frame") <- sf
  list(Y = Y, event_model = emod, X_list = X_list,
       h_true = h_true, beta_true = beta_true, sframe = sf)
}


test_that("hrfals works across HRF bases", {
  bases <- list(HRF_SPMG3, gen_hrf(hrf_bspline, N=4))
  for (b in bases) {
    dat <- simulate_cfals_wrapper_data(b)
    design <- create_cfals_design(dat$Y, dat$event_model, b)
    fit <- hrfals(dat$Y, dat$event_model, b,
                             lam_beta = 0.1, lam_h = 0.1)
    expect_equal(dim(fit$h_coeffs), c(nbasis(b), ncol(dat$Y)))
    expect_equal(dim(fit$beta_amps), c(length(dat$X_list), ncol(dat$Y)))
    recon <- design$Phi_recon_matrix %*% fit$h_coeffs
    expect_true(all(is.finite(recon)))
  }
})


test_that("cfals handles low-signal data", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3, noise_sd = 0.5, signal_scale = 0.01)
  fit <- hrfals(dat$Y, dat$event_model, HRF_SPMG3)
  expect_lt(mean(fit$gof_per_voxel), 0.2)
})

simple_cfals_data_noise <- function() {
  set.seed(123)
  n <- 50
  d <- 3
  k <- 2
  v <- 3
  h_true <- matrix(rnorm(d * v), d, v)
  beta_true <- matrix(rnorm(k * v), k, v)
  X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
  Xbig <- do.call(cbind, X_list)
  Y <- Xbig %*% as.vector(matrix(h_true, d, v) %*% t(beta_true))
  Y <- matrix(Y, n, v) + matrix(rnorm(n * v, sd = 0.01), n, v)
  phi <- diag(d)
  href <- rep(1, nrow(phi))
  list(X_list = X_list, Y = Y, Xbig = Xbig, Phi = phi, href = href)
}

test_that("cf_als_engine predictions match canonical GLM", {
  dat <- simple_cfals_data_noise()
  res <- cf_als_engine(dat$X_list, dat$Y,
                       lambda_b = 0,
                       lambda_h = 0,
                       R_mat_eff = NULL,
                       fullXtX_flag = FALSE,
                       precompute_xty_flag = TRUE,
                       Phi_recon_matrix = dat$Phi,
                       h_ref_shape_canonical = dat$href,
                       max_alt = 1)
  n <- nrow(dat$Y)
  v <- ncol(dat$Y)
  pred_cfals <- matrix(0, n, v)
  for (c in seq_along(dat$X_list)) {
    pred_cfals <- pred_cfals + (dat$X_list[[c]] %*% res$h) *
      matrix(rep(res$beta[c, ], each = n), n, v)
  }
  gamma_hat <- chol2inv(chol(crossprod(dat$Xbig))) %*% crossprod(dat$Xbig, dat$Y)
  pred_glm <- dat$Xbig %*% gamma_hat
  expect_equal(pred_cfals, pred_glm, tolerance = 1.0)
})


test_that("hrfals predictions match canonical GLM", {
  set.seed(123)
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  fit <- hrfals(dat$Y, dat$event_model, HRF_SPMG3,
                lam_beta = 0,
                lam_h = 0,
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


