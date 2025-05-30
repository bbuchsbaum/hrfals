context("simulate_simple_dataset integration")

library(fmrireg)

# generate dataset via helper and fit models across wrappers

test_that("cfals wrappers work with simulated dataset", {
  dat <- simulate_simple_data(ncond = 2, snr = 10)

  # determine fmri data matrix and sampling frame
  Y <- dat$noisy  # Use noisy data directly (no time column to remove)
  sframe <- dat$sampling_frame
  attr(Y, "sampling_frame") <- sframe

  events <- data.frame(
    onset = dat$onsets,
    condition = factor(dat$conditions),
    block = 1
  )

  emod <- event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sframe)

  nb <- nbasis(fmrireg::HRF_SPMG2)
  nvox <- ncol(Y)
  k <- length(unique(events$condition))

  fit2 <- hrfals(Y, emod, fmrireg::HRF_SPMG2,
                 lam_beta = 0, lam_h = 0)
  fit3 <- estimate_hrf_cfals(Y, emod, "hrf(condition)",
                             fmrireg::HRF_SPMG2,
                             lambda_b = 0, lambda_h = 0)

  expect_equal(dim(fit2$h_coeffs), c(nb, nvox))
  expect_true(mean(fit2$gof_per_voxel) > 0)
  expect_equal(dim(fit3$h_coeffs), c(nb, nvox))
  expect_true(mean(fit3$gof_per_voxel) > 0)
})
