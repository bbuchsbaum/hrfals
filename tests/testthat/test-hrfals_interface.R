context("hrfals interface")

library(fmrireg)

test_that("hrfals forwards arguments to hrfals wrapper", {
  dat <- simulate_cfals_wrapper_data(fmrihrf::HRF_SPMG3)
  design <- create_cfals_design(dat$Y, dat$event_model, fmrihrf::HRF_SPMG3)

  ctrl <- list(lambda_init = 0, lambda_b = 0.1, lambda_h = 0.1,
               fullXtX = TRUE, max_alt = 2)

  # Test the hrfals_from_design interface function (design-based)
  fit1 <- hrfals_from_design(dat$Y, design, method = "cf_als", control = ctrl)
  
  # Test the direct hrfals function (event_model-based)
  fit2 <- hrfals(dat$Y, design$event_model, design$hrf_basis,
                 lam_beta = 0.1, lam_h = 0.1, fullXtX = TRUE, max_alt = 2)

  expect_equal(fit1$h_coeffs, fit2$h_coeffs)
  expect_equal(fit1$beta_amps, fit2$beta_amps)
})
