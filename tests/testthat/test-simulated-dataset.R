context("simulate_simple_dataset integration")

library(fmrireg)

# generate dataset via helper and fit models across wrappers

test_that("cfals wrappers work with simulated dataset", {
  dat <- simulate_simple_data(ncond = 2, snr = 1)

  # determine fmri data matrix and sampling frame
  Y <- if (is.matrix(dat$dataset)) dat$dataset else fmrireg::get_data_matrix(dat$dataset)
  sframe <- if (inherits(dat$dataset, "fmri_dataset")) dat$dataset$sampling_frame else attr(dat$dataset, "sampling_frame")

  events <- data.frame(
    onset = dat$onsets,
    condition = factor(dat$conditions)
  )

  emod <- event_model(onset ~ hrf(condition), data = events,
                      sampling_frame = sframe)

  nb <- nbasis(fmrireg::HRF_SPMG3)
  nvox <- ncol(Y)
  k <- length(unique(events$condition))

  fit1 <- fmrireg_cfals(Y, emod, fmrireg::HRF_SPMG3,
                        lambda_b = 0.1, lambda_h = 0.1)
  fit2 <- fmrireg_hrf_cfals(Y, emod, fmrireg::HRF_SPMG3,
                            lam_beta = 0.1, lam_h = 0.1)
  fit3 <- estimate_hrf_cfals(Y, emod, "hrf(condition)",
                             fmrireg::HRF_SPMG3,
                             lambda_b = 0.1, lambda_h = 0.1)

  expect_equal(dim(fit1$h_coeffs), c(nb, nvox))
  expect_equal(dim(fit1$beta_amps), c(k, nvox))
  expect_true(mean(fit1$gof_per_voxel) > 0)
  expect_equal(dim(fit2$h_coeffs), c(nb, nvox))
  expect_true(mean(fit2$gof_per_voxel) > 0)
  expect_equal(dim(fit3$h_coeffs), c(nb, nvox))
  expect_true(mean(fit3$gof_per_voxel) > 0)
})
