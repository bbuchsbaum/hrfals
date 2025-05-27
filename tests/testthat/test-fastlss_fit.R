context("fastlss_fit object")

make_cfals_fit <- function() {
  h <- matrix(rnorm(2), 1, 2)
  beta <- matrix(rnorm(2), 1, 2)
  phi <- matrix(1, 1, 1)
  dinfo <- list(d = 1, k = 1, n = 1, v = 2)
  hrfals_fit(h, beta, "cf_als", c(beta = 1, h = 1),
             call("hrfals_fit"), fmrireg::HRF_SPMG1, "term", phi, dinfo,
             matrix(0, 1, 1))
}

test_that("fastlss_fit constructor works", {
  betas <- matrix(rnorm(6), 3, 2)
  cf_fit <- make_cfals_fit()
  fit <- fastlss_fit(betas, mode = "shared", cfals_fit = cf_fit,
                     events = data.frame(onset = 1, condition = "A"),
                     hrf_basis = fmrireg::HRF_SPMG1,
                     call = quote(test_call()))
  expect_s3_class(fit, "fastlss_fit")
  expect_equal(fit$betas, betas)
  expect_equal(fit$mode, "shared")
})


test_that("fastlss_fit print and summary work", {
  betas <- matrix(rnorm(6), 3, 2)
  cf_fit <- make_cfals_fit()
  fit <- fastlss_fit(betas, mode = "shared", cfals_fit = cf_fit,
                     events = data.frame(onset = 1, condition = "A"),
                     hrf_basis = fmrireg::HRF_SPMG1,
                     call = quote(test_call()))
  expect_output(print(fit), "fastLSS beta-series")
  expect_output(print(summary(fit)), "Summary of fastLSS fit")
})
