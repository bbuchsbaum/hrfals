library(testthat)
library(fmridesign)

context("Simple HRF Recovery Test")

test_that("ls_svd_engine recovers correct HRF shape in simple case", {
  # This is a very basic test of HRF reconstruction using ls_svd_engine
  # We generate simple synthetic data with a known HRF shape and verify recovery
  
  set.seed(123)
  
  # Simple setup: single condition, single voxel
  TR <- 2.0
  n_timepoints <- 50
  timegrid <- seq(0, (n_timepoints - 1) * TR, by = TR)
  
  # Create two different HRF shapes to test
  hrf1 <- fmrihrf::block_hrf(fmrihrf::lag_hrf(fmrihrf::HRF_GAUSSIAN, lag = 2), 3)
  hrf2 <- fmrihrf::block_hrf(fmrihrf::lag_hrf(fmrihrf::HRF_GAUSSIAN, lag = 4), 5)
  
  # Evaluate both HRFs on our time grid
  hrf1_vals <- fmrihrf::evaluate(hrf1, timegrid)
  hrf2_vals <- fmrihrf::evaluate(hrf2, timegrid)
  if (is.matrix(hrf1_vals)) hrf1_vals <- hrf1_vals[, 1]
  if (is.matrix(hrf2_vals)) hrf2_vals <- hrf2_vals[, 1]
  
  # Create simple neural signal (few impulses)
  neural_signal <- rep(0, n_timepoints)
  neural_signal[c(5, 15, 25)] <- 1  # Three events
  
  # Generate BOLD signals by convolving with each HRF
  Y1_clean <- stats::convolve(neural_signal, rev(hrf1_vals), type = "open")[1:n_timepoints]
  Y2_clean <- stats::convolve(neural_signal, rev(hrf2_vals), type = "open")[1:n_timepoints]
  
  # Add small amount of noise
  noise_level <- 0.05
  Y1_noisy <- Y1_clean + rnorm(n_timepoints, 0, noise_level * sd(Y1_clean))
  Y2_noisy <- Y2_clean + rnorm(n_timepoints, 0, noise_level * sd(Y2_clean))
  
  # Create design matrix using FIR basis
  hrf_basis <- fmrihrf::HRF_FIR
  d <- fmrihrf::nbasis(hrf_basis)
  
  # Build design matrix for single condition
  X_design <- matrix(0, n_timepoints, d)
  basis_vals <- fmrihrf::evaluate(hrf_basis, timegrid)
  
  for (j in seq_len(d)) {
    if (is.matrix(basis_vals)) {
      basis_j <- basis_vals[, j]
    } else {
      basis_j <- basis_vals
    }
    X_design[, j] <- stats::convolve(neural_signal, rev(basis_j), type = "open")[1:n_timepoints]
  }
  
  # Create reconstruction matrix and canonical reference
  Phi_recon <- reconstruction_matrix(hrf_basis, timegrid)
  h_ref_canonical <- fmridesign::evaluate(fmrihrf::HRF_SPMG1, timegrid)
  if (is.matrix(h_ref_canonical)) h_ref_canonical <- h_ref_canonical[, 1]
  h_ref_canonical <- h_ref_canonical / max(abs(h_ref_canonical))
  
  # Test 1: Run ls_svd_engine on data generated with hrf1
  result1 <- hrfals:::ls_svd_engine(
    X_list_proj = list(cond1 = X_design),
    Y_proj = matrix(Y1_noisy, ncol = 1),
    lambda_init = 1,
    Phi_recon_matrix = Phi_recon,
    h_ref_shape_canonical = h_ref_canonical
  )
  
  # Test 2: Run ls_svd_engine on data generated with hrf2
  result2 <- hrfals:::ls_svd_engine(
    X_list_proj = list(cond1 = X_design),
    Y_proj = matrix(Y2_noisy, ncol = 1),
    lambda_init = 1,
    Phi_recon_matrix = Phi_recon,
    h_ref_shape_canonical = h_ref_canonical
  )
  
  # Reconstruct estimated HRF shapes
  estimated_hrf1 <- Phi_recon %*% result1$h[, 1]
  estimated_hrf2 <- Phi_recon %*% result2$h[, 1]
  
  # Test that each estimated HRF correlates better with its true shape
  cor1_with_true1 <- cor(estimated_hrf1, hrf1_vals)
  cor1_with_true2 <- cor(estimated_hrf1, hrf2_vals)
  cor2_with_true1 <- cor(estimated_hrf2, hrf1_vals)
  cor2_with_true2 <- cor(estimated_hrf2, hrf2_vals)
  
  # The key test: each estimated HRF should correlate better with its true shape
  expect_gt(cor1_with_true1, cor1_with_true2)  # hrf1 data -> better match with hrf1
  expect_gt(cor2_with_true2, cor2_with_true1)  # hrf2 data -> better match with hrf2
  
  # Both should have reasonable correlations with their true shapes
  expect_gt(cor1_with_true1, 0.6)
  expect_gt(cor2_with_true2, 0.6)
  
  # Basic sanity checks
  expect_equal(dim(result1$h), c(d, 1))
  expect_equal(dim(result1$beta), c(1, 1))
  expect_equal(dim(result2$h), c(d, 1))
  expect_equal(dim(result2$beta), c(1, 1))
  
  # Beta coefficients should be positive (we used positive neural signals)
  expect_gt(result1$beta[1, 1], 0)
  expect_gt(result2$beta[1, 1], 0)
}) 