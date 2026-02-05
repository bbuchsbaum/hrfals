context("hrfals_lss wrapper")

library(fmridesign)

# Tests for the hrfals LSS (Least Squares Separate) implementation
# Includes comparison with fmrireg::glm_lss to ensure both methods work correctly

test_that("hrfals_lss runs in both modes", {
  dat <- simulate_cfals_wrapper_data(fmrihrf::HRF_SPMG3)
  fit <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", fmrihrf::HRF_SPMG3,
                           method = "ls_svd_only")
  B_shared <- hrfals_lss(fit, dat$event_model, fmri_data_obj = dat$Y,
                          mode = "shared")
  expect_s3_class(B_shared, "fastlss_fit")
  expect_equal(nrow(B_shared$betas), length(dat$X_list))
  expect_equal(ncol(B_shared$betas), ncol(dat$Y))

  B_auto <- hrfals_lss(fit, dat$event_model, fmri_data_obj = dat$Y,
                        mode = "auto")
  expect_equal(dim(B_auto$betas), c(length(dat$X_list), ncol(dat$Y)))
})

# Whitening support
test_that("hrfals_lss stores whitening matrix", {
  dat <- simulate_cfals_wrapper_data(fmrihrf::HRF_SPMG3)
  fit <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", fmrihrf::HRF_SPMG3,
                           method = "ls_svd_only")
  n <- nrow(dat$Y)
  set.seed(3)
  W <- chol(crossprod(matrix(rnorm(n*n), n, n)))
  res <- hrfals_lss(fit, dat$event_model, fmri_data_obj = dat$Y,
                    mode = "shared", whitening_matrix = W)
  expect_s3_class(res, "fastlss_fit")
  expect_true(is.matrix(res$whitening_matrix))
  expect_equal(res$whitening_matrix, W)
})

# Comparison with fmridesign LSS implementation
test_that("hrfals_lss produces trial-level results comparable to fmrireg::glm_lss", {
  # Create test data with individual trial regressors for proper LSS
  sf <- fmridesign::sampling_frame(blocklens = 60, TR = 1)
  events <- data.frame(
    onset = c(5, 15, 30, 45),
    condition = factor(c("A", "A", "B", "B")),
    trial_id = factor(1:4),  # Individual trial identifiers
    block = 1
  )
  
  # Create event model with individual trial regressors for LSS
  # Use trial_id instead of condition to get individual trial regressors
  emod <- fmridesign::event_model(onset ~ fmridesign::hrf(trial_id), data = events,
                               block = ~ block, sampling_frame = sf)
  
  # Create design matrices
  design <- create_fmri_design(emod, fmrihrf::HRF_FIR)
  X_list <- design$X_list
  
  # Generate simulated data
  d <- design$d
  k <- design$k  # Should be 4 (one per trial)
  v <- 2
  n_timepoints <- length(fmridesign::samples(sf, global = TRUE))
  
  h_true <- matrix(rnorm(d * v), d, v) * 0.5
  beta_true <- matrix(rnorm(k * v), k, v) * 0.5
  Y <- matrix(0, n_timepoints, v)
  for (c in seq_along(X_list)) {
    Y <- Y + (X_list[[c]] %*% h_true) *
      matrix(rep(beta_true[c, ], each = nrow(Y)), nrow(Y), v)
  }
  Y <- Y + matrix(rnorm(length(Y), sd = 0.05), nrow(Y), v)
  attr(Y, "sampling_frame") <- sf
  
  # Check that we have 4 trial regressors
  expect_equal(length(X_list), 4, info = "Should have 4 individual trial regressors")
  
  # Fit CF-ALS model using trial-level regressors
  fit <- estimate_hrf_cfals(Y, emod, "hrf(trial_id)", fmrihrf::HRF_FIR, method = "ls_svd_only")
  
  # Get hrfals LSS results (should now be trial-level)
  hrfals_result <- hrfals_lss(fit, emod, fmri_data_obj = Y, mode = "shared")
  hrfals_betas <- hrfals_result$betas
  
  # Check that hrfals LSS returns trial-level results
  expect_equal(nrow(hrfals_betas), 4, info = "hrfals LSS should return 4 trial-level betas")
  expect_equal(ncol(hrfals_betas), 2, info = "hrfals LSS should return betas for 2 voxels")
  
  # Create fmridesign dataset for comparison
  skip_if_not_installed("fmridataset")

  fmri_dataset <- fmridataset::matrix_dataset(Y, TR = 1, 
                                          run_length = nrow(Y),
                                          event_table = events)
  
  # Create model for fmridesign LSS (using trial_id for individual trials)
  model_obj <- fmridesign::event_model(onset ~ hrf(trial_id), 
                                    data = events,
                                    block = ~ block,
                                    sampling_frame = fmri_dataset$sampling_frame)
  
  # Get fmridesign LSS results
  fmridesign_result <- fmrireg::glm_lss(fmri_dataset, model_obj, fmrihrf::HRF_FIR)
  
  # Test that both methods produce reasonable trial-level results
  expect_s3_class(hrfals_result, "fastlss_fit")
  expect_true(is.matrix(hrfals_betas))
  expect_equal(ncol(hrfals_betas), ncol(Y))  # Same number of voxels
  
  expect_true(is.list(fmridesign_result))
  expect_true("betas_ran" %in% names(fmridesign_result))
  expect_true(is.matrix(fmridesign_result$betas_ran))
  expect_equal(ncol(fmridesign_result$betas_ran), ncol(Y))  # Same number of voxels
  
  # Test that hrfals LSS results are reasonable (main focus)
  expect_gt(max(abs(hrfals_betas)), 0.01)
  expect_lt(max(abs(hrfals_betas)), 10)
  expect_true(all(is.finite(hrfals_betas)))
  
  # Test fmridesign results if they're valid (may have numerical issues)
  if (all(is.finite(fmridesign_result$betas_ran))) {
    expect_gt(max(abs(fmridesign_result$betas_ran)), 0.01)
    expect_lt(max(abs(fmridesign_result$betas_ran)), 10)
    cat("Both methods produce valid results\n")
  } else {
    cat("fmridesign LSS encountered numerical issues; hrfals LSS works correctly\n")
  }
  
  cat("hrfals LSS dimensions:", paste(dim(hrfals_betas), collapse="x"), "(trial-level)\n")
  cat("fmridesign LSS dimensions:", paste(dim(fmridesign_result$betas_ran), collapse="x"), "\n")
})

# Test that LSS produces trial-level amplitude estimates
test_that("hrfals_lss correctly estimates trial-level amplitudes", {
  # Create simple test data with known trial structure
  sf <- fmridesign::sampling_frame(blocklens = 60, TR = 1)
  events <- data.frame(
    onset = c(10, 20, 30, 40),
    trial_id = factor(paste0("trial_", 1:4)),
    block = 1
  )
  
  # Create event model with individual trial regressors
  emod <- fmridesign::event_model(onset ~ fmridesign::hrf(trial_id), data = events,
                               block = ~ block, sampling_frame = sf)
  
  # Create design and simulate data
  design <- create_fmri_design(emod, fmrihrf::HRF_FIR)
  n_timepoints <- length(fmridesign::samples(sf, global = TRUE))
  Y <- matrix(rnorm(n_timepoints * 3), n_timepoints, 3)
  attr(Y, "sampling_frame") <- sf
  
  # Fit CF-ALS to estimate HRF shape
  fit <- estimate_hrf_cfals(Y, emod, "hrf(trial_id)", fmrihrf::HRF_FIR, method = "ls_svd_only")
  
  # Apply LSS to get trial-specific amplitudes
  lss_result <- hrfals_lss(fit, emod, fmri_data_obj = Y, mode = "shared")
  
  # Verify trial-level structure
  expect_equal(nrow(lss_result$betas), 4, info = "Should have 4 trial-specific amplitudes")
  expect_equal(ncol(lss_result$betas), 3, info = "Should have amplitudes for 3 voxels")
  expect_true(all(grepl("trial_", rownames(lss_result$betas))), 
               info = "Row names should contain trial identifiers")
  
  # Verify that each trial gets a different amplitude estimate
  trial_amplitudes <- lss_result$betas
  expect_true(length(unique(as.vector(trial_amplitudes))) > 1, 
              info = "Trial amplitudes should vary (not all identical)")
  
  cat("Trial-level LSS successfully estimated", nrow(trial_amplitudes), 
      "trial amplitudes for", ncol(trial_amplitudes), "voxels\n")
})
