context("hrfals_fit class")

test_that("constructor creates object with correct class", {
  h <- matrix(0, 2, 1)
  b <- matrix(0, 1, 1)
  dinfo <- list(d = 2, k = 1, n = 1, v = 1, fullXtX = FALSE)
  phi <- diag(2)
  fit <- hrfals_fit(h, b, "test", c(init = 0), call("foo"),
                    fmrireg::HRF_SPMG1, "term", phi, dinfo, matrix(0,1,1))
  expect_s3_class(fit, "hrfals_fit")
  expect_equal(fit$method_used, "test")
  expect_equal(fit$target_event_term_name, "term")
  expect_equal(fit$phi_recon_matrix, phi)
})

test_that("print and summary methods work", {
  h <- matrix(0, 2, 1)
  b <- matrix(0, 1, 1)
  dinfo <- list(d = 2, k = 1, n = 1, v = 1, fullXtX = FALSE)
  phi <- diag(2)
  fit <- hrfals_fit(h, b, "test", c(init = 0), call("foo"),
                    fmrireg::HRF_SPMG1, "term", phi, dinfo, matrix(0,1,1))
  expect_output(print(fit), "hrfals Fit")
  s <- summary(fit)
  expect_s3_class(s, "summary.hrfals_fit")
  expect_equal(s$lambdas, c(init = 0))
  expect_invisible(plot(fit, vox = 1))
})
