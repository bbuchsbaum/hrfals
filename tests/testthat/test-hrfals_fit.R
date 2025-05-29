context("hrfals_fit object")

test_that("hrfals_fit constructor works", {
  h <- matrix(rnorm(6), 3, 2)
  beta <- matrix(rnorm(4), 2, 2)
  phi <- matrix(rnorm(15), 5, 3)
  dinfo <- list(d = 3, k = 2, n = 5, v = 2,
                predictor_means = c(0, 0), predictor_sds = c(1, 1))
  fit <- hrfals_fit(h, beta, "cf_als", c(beta = 1, h = 1),
                    quote(test_call()), fmrireg::HRF_SPMG3, "term", phi, dinfo, matrix(0,1,1),
                    beta_penalty = list(l1 = 0.1, alpha = 1, warm_start = TRUE))
  expect_s3_class(fit, "hrfals_fit")
  expect_equal(fit$h_coeffs, h)
  expect_equal(fit$beta_amps, beta)
  expect_equal(fit$lambdas$beta_penalty$l1, 0.1)
  expect_equal(fit$design_info$predictor_sds, c(1, 1))
})

test_that("print.hrfals_fit doesn't error", {
  h <- matrix(rnorm(6), 3, 2)
  beta <- matrix(rnorm(4), 2, 2)
  phi <- matrix(rnorm(15), 5, 3)
  dinfo <- list(d = 3, k = 2, n = 5, v = 2,
                predictor_means = c(0, 0), predictor_sds = c(1, 1))
  fit <- hrfals_fit(h, beta, "cf_als", c(beta = 1, h = 1),
                    quote(test_call()), fmrireg::HRF_SPMG3, "term", phi, dinfo, matrix(0,1,1),
                    beta_penalty = list(l1 = 0.1, alpha = 1, warm_start = TRUE))
  expect_output(print(fit), "hrfals Fit")
})
