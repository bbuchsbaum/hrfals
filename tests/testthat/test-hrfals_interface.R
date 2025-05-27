context("hrfals interface")

library(fmrireg)

test_that("hrfals forwards arguments to hrfals wrapper", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  design <- create_cfals_design(dat$Y, dat$event_model, HRF_SPMG3)

  ctrl <- list(lambda_init = 0, lambda_b = 0.1, lambda_h = 0.1,
               fullXtX = TRUE, max_alt = 2)

  fit1 <- hrfals(dat$Y, design, method = "cf_als", control = ctrl)
  fit2 <- hrfals(dat$Y, design$event_model, design$hrf_basis,
                 lam_beta = 0.1, lam_h = 0.1,
                 fullXtX = TRUE, max_alt = 2)

  expect_equal(fit1$h_coeffs, fit2$h_coeffs)
  expect_equal(fit1$beta_amps, fit2$beta_amps)
})
