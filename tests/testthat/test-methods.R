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

make_fmrireg_fit <- function() {
  h <- matrix(rnorm(4), 2, 2)
  beta <- matrix(rnorm(2), 1, 2)
  recon <- matrix(rnorm(4), 2, 2)
  dinfo <- list(d = 2, k = 1, n = 2, v = 2)
  fmrireg_cfals_fit(h, beta, "cf_als", c(beta = 1, h = 1),
                    call("fmrireg_cfals"), fmrireg::HRF_SPMG1,
                    dinfo, matrix(0, 1, 1), recon)
}


test_that("hrfals_fit methods run", {
  fit <- make_hrfals_fit()
  expect_output(print(fit), "hrfals Fit")
  expect_output(print(summary(fit)), "Summary of hrfals Fit")
  pdf(NULL)
  expect_silent(plot(fit))
  dev.off()
})


test_that("fmrireg_cfals_fit methods run", {
  fit <- make_fmrireg_fit()
  expect_output(print(fit), "fmrireg CF-ALS Fit")
  expect_output(print(summary(fit)), "Summary of fmrireg CF-ALS Fit")
  pdf(NULL)
  expect_silent(plot(fit))
  dev.off()
})

