library(testthat)
library(fmridesign)

context("HRF Shape and Amplitude Recovery")

test_that("ls_svd_engine recovers HRF shape and amplitude with same shape", {
  # Test that ls_svd_engine can recover both the correct HRF shape AND amplitude
  # when we have two conditions with the SAME HRF shape but different amplitudes
  # This is what ls_svd_engine is designed for: rank-1 decomposition with shared HRF
  
  set.seed(456)
  
  # Experimental setup
  TR <- 2.0
  n_timepoints <- 80
  timegrid <- seq(0, (n_timepoints - 1) * TR, by = TR)
  
  # Same HRF shape for both conditions (this is key for ls_svd_engine)
  hrf_shape <- fmrihrf::HRF_SPMG1
  
  # Define different amplitudes for each condition
  amplitude1 <- 2.5   # Condition 1: high amplitude
  amplitude2 <- 1.0   # Condition 2: lower amplitude
  
  # Create event onsets for each condition
  onsets_cond1 <- c(10, 40, 70, 100)   # Condition 1 events (in seconds)
  onsets_cond2 <- c(20, 50, 80, 110)   # Condition 2 events (in seconds)
  
  # Create regressors with same HRF shape but different amplitudes
  reg1 <- fmridesign::regressor(onsets = onsets_cond1, hrf = hrf_shape, 
                            amplitude = amplitude1, duration = 0)
  reg2 <- fmridesign::regressor(onsets = onsets_cond2, hrf = hrf_shape, 
                            amplitude = amplitude2, duration = 0)
  
  # Generate clean BOLD signals using evaluate (this handles convolution automatically)
  Y1_clean <- fmridesign::evaluate(reg1, timegrid)
  Y2_clean <- fmridesign::evaluate(reg2, timegrid)
  if (is.matrix(Y1_clean)) Y1_clean <- Y1_clean[, 1]
  if (is.matrix(Y2_clean)) Y2_clean <- Y2_clean[, 1]
  
  # Combine signals from both conditions (single voxel responds to both)
  Y_combined_clean <- Y1_clean + Y2_clean
  
  # Add small amount of noise
  noise_level <- 0.05
  Y_noisy <- Y_combined_clean + rnorm(n_timepoints, 0, noise_level * sd(Y_combined_clean))
  
  # Create design matrices using FIR basis
  hrf_basis <- fmrihrf::HRF_FIR
  d <- fmrihrf::nbasis(hrf_basis)
  
  # Create neural signals (stick functions) for each condition
  neural_signal1 <- rep(0, n_timepoints)
  neural_signal2 <- rep(0, n_timepoints)
  
  for (onset in onsets_cond1) {
    idx <- which.min(abs(timegrid - onset))
    if (idx <= length(neural_signal1)) neural_signal1[idx] <- 1
  }
  
  for (onset in onsets_cond2) {
    idx <- which.min(abs(timegrid - onset))
    if (idx <= length(neural_signal2)) neural_signal2[idx] <- 1
  }
  
  # Build design matrices for both conditions
  X_design1 <- matrix(0, n_timepoints, d)
  X_design2 <- matrix(0, n_timepoints, d)
  basis_vals <- fmridesign::evaluate(hrf_basis, timegrid)
  
  for (j in seq_len(d)) {
    if (is.matrix(basis_vals)) {
      basis_j <- basis_vals[, j]
    } else {
      basis_j <- basis_vals
    }
    X_design1[, j] <- stats::convolve(neural_signal1, rev(basis_j), type = "open")[1:n_timepoints]
    X_design2[, j] <- stats::convolve(neural_signal2, rev(basis_j), type = "open")[1:n_timepoints]
  }
  
  # Create reconstruction matrix and canonical reference
  Phi_recon <- reconstruction_matrix(hrf_basis, timegrid)
  h_ref_canonical <- fmridesign::evaluate(fmrihrf::HRF_SPMG1, timegrid)
  if (is.matrix(h_ref_canonical)) h_ref_canonical <- h_ref_canonical[, 1]
  h_ref_canonical <- h_ref_canonical / max(abs(h_ref_canonical))
  
  # Run both ls_svd_engine and ls_svd_1als_engine for comparison
  result_svd <- ls_svd_engine(
    X_list_proj = list(cond1 = X_design1, cond2 = X_design2),
    Y_proj = matrix(Y_noisy, ncol = 1),
    lambda_init = 1,
    Phi_recon_matrix = Phi_recon,
    h_ref_shape_canonical = h_ref_canonical,
    epsilon_svd = 1e-8,
    epsilon_scale = 1e-8
  )
  
  result_als <- ls_svd_1als_engine(
    X_list_proj = list(cond1 = X_design1, cond2 = X_design2),
    Y_proj = matrix(Y_noisy, ncol = 1),
    lambda_init = 1,
    lambda_b = 10,
    lambda_h = 1,
    Phi_recon_matrix = Phi_recon,
    h_ref_shape_canonical = h_ref_canonical,
    epsilon_svd = 1e-8,
    epsilon_scale = 1e-8
  )
  
  # Test 1: Shape Recovery
  # Reconstruct estimated HRF shapes for both methods
  estimated_hrf_svd <- Phi_recon %*% result_svd$h[, 1]
  estimated_hrf_als <- Phi_recon %*% result_als$h[, 1]
  true_hrf <- fmridesign::evaluate(hrf_shape, timegrid)
  if (is.matrix(true_hrf)) true_hrf <- true_hrf[, 1]
  
  # Normalize for fair comparison
  true_hrf_norm <- true_hrf / max(abs(true_hrf))
  estimated_hrf_svd_norm <- estimated_hrf_svd / max(abs(estimated_hrf_svd))
  estimated_hrf_als_norm <- estimated_hrf_als / max(abs(estimated_hrf_als))
  
  # Both should have high correlation with true HRF shape
  correlation_svd <- cor(estimated_hrf_svd_norm, true_hrf_norm)
  correlation_als <- cor(estimated_hrf_als_norm, true_hrf_norm)
  
  expect_gt(correlation_svd, 0.8)
  expect_gt(correlation_als, 0.8)
  
  # ALS should do at least as well as SVD (often better)
  expect_gte(correlation_als, correlation_svd - 0.05)  # Allow small tolerance
  
  # Test 2: Amplitude Recovery
  # Compare amplitude recovery for both methods
  beta1_svd <- result_svd$beta[1, 1]  # SVD: Condition 1 amplitude estimate
  beta2_svd <- result_svd$beta[2, 1]  # SVD: Condition 2 amplitude estimate
  beta1_als <- result_als$beta[1, 1]  # ALS: Condition 1 amplitude estimate
  beta2_als <- result_als$beta[2, 1]  # ALS: Condition 2 amplitude estimate
  
  # All betas should be positive (we used positive amplitudes)
  expect_gt(beta1_svd, 0)
  expect_gt(beta2_svd, 0)
  expect_gt(beta1_als, 0)
  expect_gt(beta2_als, 0)
  
  # The higher amplitude condition should have higher beta for both methods
  expect_gt(beta1_svd, beta2_svd)
  expect_gt(beta1_als, beta2_als)
  
  # Compare amplitude ratio recovery
  estimated_ratio_svd <- beta1_svd / beta2_svd
  estimated_ratio_als <- beta1_als / beta2_als
  true_ratio <- amplitude1 / amplitude2  # 2.5 / 1.0 = 2.5
  
  # Both should be reasonably accurate for same HRF shape
  ratio_error_svd <- abs(estimated_ratio_svd - true_ratio) / true_ratio
  ratio_error_als <- abs(estimated_ratio_als - true_ratio) / true_ratio
  
  expect_lt(ratio_error_svd, 0.8)  # Within 80% of true ratio
  expect_lt(ratio_error_als, 0.8)  # Within 80% of true ratio
  
  # ALS should do at least as well as SVD for amplitude recovery
  expect_lte(ratio_error_als, ratio_error_svd + 0.1)  # Allow small tolerance
  
  # Test 3: Basic sanity checks for both methods
  expect_equal(dim(result_svd$h), c(d, 1))
  expect_equal(dim(result_svd$beta), c(2, 1))  # Two conditions
  expect_equal(nrow(result_svd$Gamma_hat), d * 2)  # d coefficients × 2 conditions
  
  expect_equal(dim(result_als$h), c(d, 1))
  expect_equal(dim(result_als$beta), c(2, 1))  # Two conditions
  expect_equal(nrow(result_als$Gamma_hat), d * 2)  # d coefficients × 2 conditions
})

test_that("ls_svd_engine handles amplitude differences with same HRF shape", {
  # Simpler test: same HRF shape, different amplitudes
  # This isolates amplitude recovery from shape recovery
  
  set.seed(789)
  
  TR <- 2.0
  n_timepoints <- 60
  timegrid <- seq(0, (n_timepoints - 1) * TR, by = TR)
  
  # Same HRF shape for both conditions
  hrf_shape <- fmrihrf::HRF_SPMG1
  
  # Different amplitudes
  amplitude1 <- 3.0
  amplitude2 <- 1.0
  
  # Different event timings
  onsets_cond1 <- c(10, 40, 70)
  onsets_cond2 <- c(20, 50, 80)
  
  # Create regressors
  reg1 <- fmridesign::regressor(onsets = onsets_cond1, hrf = hrf_shape, 
                            amplitude = amplitude1, duration = 0)
  reg2 <- fmridesign::regressor(onsets = onsets_cond2, hrf = hrf_shape, 
                            amplitude = amplitude2, duration = 0)
  
  # Generate signals
  Y1_clean <- fmridesign::evaluate(reg1, timegrid)
  Y2_clean <- fmridesign::evaluate(reg2, timegrid)
  if (is.matrix(Y1_clean)) Y1_clean <- Y1_clean[, 1]
  if (is.matrix(Y2_clean)) Y2_clean <- Y2_clean[, 1]
  
  Y_combined <- Y1_clean + Y2_clean
  Y_noisy <- Y_combined + rnorm(n_timepoints, 0, 0.05 * sd(Y_combined))
  
  # Create design matrices (same as before)
  hrf_basis <- fmrihrf::HRF_FIR
  d <- fmrihrf::nbasis(hrf_basis)
  
  neural_signal1 <- rep(0, n_timepoints)
  neural_signal2 <- rep(0, n_timepoints)
  
  for (onset in onsets_cond1) {
    idx <- which.min(abs(timegrid - onset))
    if (idx <= length(neural_signal1)) neural_signal1[idx] <- 1
  }
  
  for (onset in onsets_cond2) {
    idx <- which.min(abs(timegrid - onset))
    if (idx <= length(neural_signal2)) neural_signal2[idx] <- 1
  }
  
  X_design1 <- matrix(0, n_timepoints, d)
  X_design2 <- matrix(0, n_timepoints, d)
  basis_vals <- fmridesign::evaluate(hrf_basis, timegrid)
  
  for (j in seq_len(d)) {
    if (is.matrix(basis_vals)) {
      basis_j <- basis_vals[, j]
    } else {
      basis_j <- basis_vals
    }
    X_design1[, j] <- stats::convolve(neural_signal1, rev(basis_j), type = "open")[1:n_timepoints]
    X_design2[, j] <- stats::convolve(neural_signal2, rev(basis_j), type = "open")[1:n_timepoints]
  }
  
  Phi_recon <- reconstruction_matrix(hrf_basis, timegrid)
  h_ref_canonical <- fmridesign::evaluate(fmrihrf::HRF_SPMG1, timegrid)
  if (is.matrix(h_ref_canonical)) h_ref_canonical <- h_ref_canonical[, 1]
  h_ref_canonical <- h_ref_canonical / max(abs(h_ref_canonical))
  
  # Run both engines for comparison
  result_svd <- ls_svd_engine(
    X_list_proj = list(cond1 = X_design1, cond2 = X_design2),
    Y_proj = matrix(Y_noisy, ncol = 1),
    lambda_init = 1,
    Phi_recon_matrix = Phi_recon,
    h_ref_shape_canonical = h_ref_canonical
  )
  
  result_als <- ls_svd_1als_engine(
    X_list_proj = list(cond1 = X_design1, cond2 = X_design2),
    Y_proj = matrix(Y_noisy, ncol = 1),
    lambda_init = 1,
    lambda_b = 10,
    lambda_h = 1,
    Phi_recon_matrix = Phi_recon,
    h_ref_shape_canonical = h_ref_canonical
  )
  
  # Test amplitude recovery for both methods
  beta1_svd <- result_svd$beta[1, 1]
  beta2_svd <- result_svd$beta[2, 1]
  beta1_als <- result_als$beta[1, 1]
  beta2_als <- result_als$beta[2, 1]
  
  expect_gt(beta1_svd, 0)
  expect_gt(beta2_svd, 0)
  expect_gt(beta1_als, 0)
  expect_gt(beta2_als, 0)
  
  # Since amplitude1 > amplitude2, beta1 should be > beta2 for both methods
  expect_gt(beta1_svd, beta2_svd)
  expect_gt(beta1_als, beta2_als)
  
  # Check amplitude ratio for both methods
  estimated_ratio_svd <- beta1_svd / beta2_svd
  estimated_ratio_als <- beta1_als / beta2_als
  true_ratio <- amplitude1 / amplitude2  # 3.0 / 1.0 = 3.0
  
  ratio_error_svd <- abs(estimated_ratio_svd - true_ratio) / true_ratio
  ratio_error_als <- abs(estimated_ratio_als - true_ratio) / true_ratio
  
  expect_lt(ratio_error_svd, 0.7)  # Within 70% for this case
  expect_lt(ratio_error_als, 0.7)  # Within 70% for this case
  
  # ALS should do at least as well as SVD
  expect_lte(ratio_error_als, ratio_error_svd + 0.1)
}) 