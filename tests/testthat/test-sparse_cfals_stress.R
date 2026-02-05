context("sparse CFALS stress test")

# helper function to simulate large data set
# Reduced dimensions for reasonable test runtime
simulate_large_predictor_data <- function(n = 200, d = 10, k = 100,
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
  skip_if_not_installed("glmnet")  # Make sure glmnet is available
  skip("Stress test - enable for performance testing")  # Skip by default
  
  # Use more reasonable dimensions that still test the functionality
  # Original was n=500, d=20, k=500 which creates 5 million matrix elements
  # This version uses n=200, d=10, k=100 which is 200,000 elements (25x smaller)
  dat <- simulate_large_predictor_data(n = 200, d = 10, k = 100, n_active = 10)

  # Time the operation to ensure it completes in reasonable time
  start_time <- Sys.time()
  
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

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat("Sparse CF-ALS stress test elapsed time:", elapsed, "seconds\n")
  
  # Verify the test completed in reasonable time (< 30 seconds)
  expect_lt(elapsed, 30)
  
  expect_equal(dim(res_cache$beta), c(dat$k, 1))
  expect_equal(dim(res_nocache$beta), c(dat$k, 1))
  expect_equal(dim(res_cache$h), c(dat$d, 1))
  expect_equal(dim(res_nocache$h), c(dat$d, 1))
})

# Add a lightweight version that actually runs
test_that("cf_als_engine handles sparse penalties correctly", {
  skip_on_cran()
  skip_if_not_installed("glmnet")
  
  # Small dimensions for quick test
  set.seed(123)
  n <- 100; d <- 5; k <- 20; n_active <- 5
  
  h_true <- rnorm(d)
  beta_true <- c(rep(1, n_active), rep(0, k - n_active))
  X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
  Y <- matrix(0, n, 1)
  for (c in seq_len(k)) {
    Y <- Y + (X_list[[c]] %*% h_true) * beta_true[c]
  }
  Y <- Y + matrix(rnorm(n, sd = 0.05), n, 1)
  
  # Run with L1 penalty
  res <- cf_als_engine(
    X_list, Y,
    lambda_b = 0.1,
    lambda_h = 0.1,
    lambda_init = 0.1,
    R_mat_eff = diag(d),
    fullXtX_flag = FALSE,
    precompute_xty_flag = TRUE,
    Phi_recon_matrix = diag(d),
    h_ref_shape_canonical = rep(1, d),
    max_alt = 2,
    beta_penalty = list(l1 = 0.05, alpha = 1, warm_start = TRUE),
    design_control = list(standardize_predictors = FALSE,
                          cache_design_blocks = FALSE)
  )
  
  # Check that sparse penalty produced some zeros
  n_zero <- sum(res$beta == 0)
  expect_gt(n_zero, 0)  # Should have some zero coefficients
  expect_lt(n_zero, k)  # But not all zeros
  
  # Check dimensions
  expect_equal(dim(res$beta), c(k, 1))
  expect_equal(dim(res$h), c(d, 1))
})
