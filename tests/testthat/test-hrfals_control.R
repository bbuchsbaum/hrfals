context("hrfals control utilities")

library(fmrireg)

# ensure default control list has expected names

test_that("hrfals_control_defaults returns expected fields", {
  defs <- hrfals_control_defaults()
  expect_type(defs, "list")
  expect_setequal(names(defs),
                  c("lambda_init", "lambda_b", "lambda_h", "lambda_joint", "R_mat",
                    "fullXtX", "precompute_xty_flag", "max_alt"))
})

# verify that overrides are respected by hrfals()

test_that("hrfals_from_design merges control overrides", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  design <- create_cfals_design(dat$Y, dat$event_model, HRF_SPMG3)
  fit <- hrfals_from_design(dat$Y, design, control = list(lambda_b = 0.5, max_alt = 2))
  expect_equal(unname(fit$lambdas["beta"]), 0.5)
  expect_equal(fit$design_info$fullXtX, FALSE)
})

# missing components in design should error

test_that("hrfals_from_design checks required design components", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  bad_design <- list(event_model = dat$event_model)
  expect_error(hrfals_from_design(dat$Y, bad_design),
               "'design' must contain 'event_model' and 'hrf_basis'")
})
