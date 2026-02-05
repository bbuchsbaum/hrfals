library(testthat)
library(fmridesign)

context("HRF Shape Recovery")

#' Generate synthetic data with known HRF shapes for testing recovery
#'
#' Creates simple synthetic BOLD data where we know the true HRF shape,
#' then tests if ls_svd_engine can recover it correctly.
generate_hrf_recovery_data <- function(true_hrf_func, noise_level = 0.1) {
  set.seed(42)
  
  # Simple experimental setup
  TR <- 2.0
  n_timepoints <- 100
  timegrid <- seq(0, (n_timepoints - 1) * TR, by = TR)
  
  # Create simple event design - single condition with a few events
  event_times <- c(20, 60, 120, 160)  # Event onsets in seconds
  
  # Create neural signal (stick function at event times)
  neural_signal <- rep(0, n_timepoints)
  for (onset in event_times) {
    idx <- which.min(abs(timegrid - onset))
    if (idx <= length(neural_signal)) {
      neural_signal[idx] <- 1
    }
  }
  
  # Evaluate true HRF on the time grid
  true_hrf_vals <- fmridesign::evaluate(true_hrf_func, timegrid)
  if (is.matrix(true_hrf_vals)) true_hrf_vals <- true_hrf_vals[, 1]
  
  # Convolve neural signal with true HRF to get clean BOLD signal
  Y_clean <- stats::convolve(neural_signal, rev(true_hrf_vals), type = "open")[1:n_timepoints]
  
  # Add a small amount of noise
  Y_noisy <- Y_clean + rnorm(n_timepoints, 0, noise_level * sd(Y_clean))
  
  # Create design matrix using FIR basis
  hrf_basis <- fmrihrf::HRF_FIR
  d <- fmrihrf::nbasis(hrf_basis)
  
  # Create design matrix for single condition
  X_design <- matrix(0, n_timepoints, d)
  
  # For each basis function, convolve neural signal
  for (j in seq_len(d)) {
    # Get j-th basis function values
    basis_vals <- fmridesign::evaluate(hrf_basis, timegrid)
    if (is.matrix(basis_vals)) {
      basis_j <- basis_vals[, j]
    } else {
      basis_j <- basis_vals
    }
    
    # Convolve neural signal with this basis function
    X_design[, j] <- stats::convolve(neural_signal, rev(basis_j), type = "open")[1:n_timepoints]
  }
  
  # Create reconstruction matrix and canonical reference
  Phi_recon <- reconstruction_matrix(hrf_basis, timegrid)
  h_ref_canonical <- fmridesign::evaluate(fmrihrf::HRF_SPMG1, timegrid)
  if (is.matrix(h_ref_canonical)) h_ref_canonical <- h_ref_canonical[, 1]
  h_ref_canonical <- h_ref_canonical / max(abs(h_ref_canonical))
  
  list(
    X_list = list(cond1 = X_design),
    Y = matrix(Y_noisy, ncol = 1),  # Single voxel
    Phi_recon = Phi_recon,
    h_ref_canonical = h_ref_canonical,
    true_hrf_vals = true_hrf_vals,
    timegrid = timegrid,
    true_hrf_func = true_hrf_func
  )
}

test_that("ls_svd_engine recovers canonical HRF shape better than alternative", {
  # Test 1: Generate data with canonical HRF (fmrihrf::HRF_SPMG1)
  data_canonical <- generate_hrf_recovery_data(fmrihrf::HRF_SPMG1, noise_level = 0.05)
  
  # Run ls_svd_engine
  result_canonical <- hrfals:::ls_svd_engine(
    X_list_proj = data_canonical$X_list,
    Y_proj = data_canonical$Y,
    lambda_init = 1,
    Phi_recon_matrix = data_canonical$Phi_recon,
    h_ref_shape_canonical = data_canonical$h_ref_canonical,
    epsilon_svd = 1e-8,
    epsilon_scale = 1e-8
  )
  
  # Reconstruct estimated HRF shape
  estimated_hrf_canonical <- data_canonical$Phi_recon %*% result_canonical$h[, 1]
  
  # Test 2: Generate data with Gaussian HRF (different shape)
  data_gaussian <- generate_hrf_recovery_data(fmrihrf::HRF_GAUSSIAN, noise_level = 0.05)
  
  # Run ls_svd_engine on Gaussian data
  result_gaussian <- hrfals:::ls_svd_engine(
    X_list_proj = data_gaussian$X_list,
    Y_proj = data_gaussian$Y,
    lambda_init = 1,
    Phi_recon_matrix = data_gaussian$Phi_recon,
    h_ref_shape_canonical = data_gaussian$h_ref_canonical,
    epsilon_svd = 1e-8,
    epsilon_scale = 1e-8
  )
  
  # Reconstruct estimated HRF shape
  estimated_hrf_gaussian <- data_gaussian$Phi_recon %*% result_gaussian$h[, 1]
  
  # Calculate correlations with true shapes
  # For canonical data, estimated HRF should correlate better with canonical than Gaussian
  cor_canonical_with_canonical <- cor(estimated_hrf_canonical, data_canonical$true_hrf_vals)
  cor_canonical_with_gaussian <- cor(estimated_hrf_canonical, data_gaussian$true_hrf_vals)
  
  # For Gaussian data, estimated HRF should correlate better with Gaussian than canonical
  cor_gaussian_with_gaussian <- cor(estimated_hrf_gaussian, data_gaussian$true_hrf_vals)
  cor_gaussian_with_canonical <- cor(estimated_hrf_gaussian, data_canonical$true_hrf_vals)
  
  # Test that the engine recovers the correct shape better than the wrong shape
  expect_gt(cor_canonical_with_canonical, cor_canonical_with_gaussian)
  
  expect_gt(cor_gaussian_with_gaussian, cor_gaussian_with_canonical)
  
  # Both correlations with correct shapes should be reasonably high
  expect_gt(cor_canonical_with_canonical, 0.7)
  
  expect_gt(cor_gaussian_with_gaussian, 0.7)
})

test_that("ls_svd_engine produces reasonable HRF estimates", {
  # Generate data with canonical HRF
  data <- generate_hrf_recovery_data(fmrihrf::HRF_SPMG1, noise_level = 0.1)
  
  # Run ls_svd_engine
  result <- hrfals:::ls_svd_engine(
    X_list_proj = data$X_list,
    Y_proj = data$Y,
    lambda_init = 1,
    Phi_recon_matrix = data$Phi_recon,
    h_ref_shape_canonical = data$h_ref_canonical,
    epsilon_svd = 1e-8,
    epsilon_scale = 1e-8
  )
  
  # Check that result has expected structure
  expect_true(is.list(result))
  expect_true("h" %in% names(result))
  expect_true("beta" %in% names(result))
  expect_true("Gamma_hat" %in% names(result))
  
  # Check dimensions
  expect_equal(nrow(result$h), fmrihrf::nbasis(fmrihrf::HRF_FIR))
  expect_equal(ncol(result$h), 1)  # Single voxel
  expect_equal(nrow(result$beta), 1)  # Single condition
  expect_equal(ncol(result$beta), 1)  # Single voxel
  
  # Reconstruct HRF and check it's reasonable
  estimated_hrf <- data$Phi_recon %*% result$h[, 1]
  
  # HRF should have positive peak (after sign alignment)
  expect_gt(max(estimated_hrf), 0)
  
  # HRF should return to baseline (approximately zero at end)
  tail_values <- tail(estimated_hrf, 3)
  expect_lt(mean(abs(tail_values)), 0.5)
  
  # Beta coefficient should be positive (we used positive neural signal)
  expect_gt(result$beta[1, 1], 0)
}) 