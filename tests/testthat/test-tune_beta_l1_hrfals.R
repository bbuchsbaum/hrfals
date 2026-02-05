context("tune_beta_l1_hrfals")

library(fmridesign)

simulate_small_data <- function() {
  sf <- sampling_frame(blocklens = 20, TR = 1)
  events <- data.frame(onset = c(5, 15), condition = factor(c("A", "B")), block = 1)
  emod <- event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  design <- create_fmri_design(emod, fmrihrf::HRF_SPMG3)
  X_list <- design$X_list
  d <- design$d
  k <- design$k
  v <- 2
  n_timepoints <- length(samples(sf, global = TRUE))
  h_true <- matrix(rnorm(d * v), d, v) * 0.5
  beta_true <- matrix(rnorm(k * v), k, v) * 0.5
  Y <- matrix(0, n_timepoints, v)
  for (c in seq_along(X_list)) {
    Y <- Y + (X_list[[c]] %*% h_true) *
      matrix(rep(beta_true[c, ], each = nrow(Y)), nrow(Y), v)
  }
  Y <- Y + matrix(rnorm(length(Y), sd = 0.05), nrow(Y), v)
  attr(Y, "sampling_frame") <- sf
  list(Y = Y, emod = emod)
}


test_that("tune_beta_l1_hrfals selects best penalty", {
  dat <- simulate_small_data()
  res <- tune_beta_l1_hrfals(dat$Y, dat$emod, fmrihrf::HRF_SPMG3,
                             l1_grid = c(0, 0.01),
                             n_outer_iterations_cfals = 1,
                             other_hrfals_args = list(lam_beta = 0.1,
                                                       lam_h = 0.1,
                                                       max_alt = 1),
                             seed = 123)
  expect_s3_class(res, "data.frame")
  expect_equal(attr(res, "best_l1"), res$l1[which.min(res$test_mse)])
  expect_true(res$best[which.min(res$test_mse)])
})
