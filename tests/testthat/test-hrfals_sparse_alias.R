context("hrfals_sparse alias")

library(fmrireg)

# simple check that hrfals_sparse forwards arguments to hrfals with defaults

test_that("hrfals_sparse replicates hrfals with sparse defaults", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  fit1 <- hrfals_sparse(dat$Y, dat$event_model, HRF_SPMG3,
                        lam_beta = 0.1, lam_h = 0.1, max_alt = 1)
  fit2 <- hrfals(dat$Y, dat$event_model, HRF_SPMG3,
                 beta_penalty = list(l1 = 0.05, alpha = 1, warm_start = TRUE),
                 design_control = list(standardize_predictors = TRUE,
                                       cache_design_blocks = TRUE),
                 lam_beta = 0.1, lam_h = 0.1, max_alt = 1)
  expect_equal(fit1$h_coeffs, fit2$h_coeffs)
  expect_equal(fit1$beta_amps, fit2$beta_amps)
})
