context("print, summary, and plot methods")

set.seed(1)

make_hrfals_fit <- function() {
  h <- matrix(rnorm(4), 2, 2)
  beta <- matrix(rnorm(2), 1, 2)
  phi <- diag(2)
  dinfo <- list(d = 2, k = 1, n = 2, v = 2)
  hrfals_fit(h, beta, "cf_als", c(beta = 1, h = 1),
             call("hrfals_fit"), fmrireg::HRF_SPMG1, "term",
             phi, dinfo, matrix(0, 1, 1))
}




test_that("hrfals_fit methods run", {
  fit <- make_hrfals_fit()
  expect_output(print(fit), "hrfals Fit")
  expect_output(print(summary(fit)), "Summary of hrfals Fit")
  pdf(NULL)
  expect_silent(plot(fit))
  dev.off()
})



test_that("tidy, glance and autoplot for hrfals_fit", {
  fit <- make_hrfals_fit()
  td <- tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), ncol(fit$h_coeffs))
  gl <- glance(fit)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
  pdf(NULL)
  expect_silent(autoplot(fit))
  dev.off()
})

