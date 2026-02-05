context("normalize_and_align_hrf and predict")

library(fmridesign)

# Test sign flipping and zeroing in normalize_and_align_hrf

test_that("normalize_and_align_hrf flips sign and zeroes small scale", {
  H <- matrix(c(-0.5, -1, -1e-8, 1e-8), 2, 2)
  B <- matrix(c(2, 3), 1, 2)
  Phi <- diag(2)
  href <- c(1, 1)
  Yp <- matrix(0, 2, 2)
  X_list <- list(cond1 = matrix(0, 2, 2))
  res <- normalize_and_align_hrf(H, B, Phi, href, epsilon_scale = 1e-6,
                                 Y_proj = Yp, X_list_proj = X_list)
  expect_equal(res$h[, 1], c(0.5, 1))
  expect_equal(unname(res$beta[1, 1]), -2)
  expect_true(all(res$h[, 2] == 0))
  expect_true(all(res$beta[, 2] == 0))
})

# Test behaviour of predict.hrfals_fit

test_that("predict.hrfals_fit returns negative residuals", {
  h <- matrix(rnorm(4), 2, 2)
  beta <- matrix(rnorm(2), 1, 2)
  phi <- diag(2)
  dinfo <- list(d = 2, k = 1, n = 1, v = 2,
                predictor_means = 0, predictor_sds = 1)
  resids <- matrix(c(0.1, -0.2), 1, 2)
  fit <- hrfals_fit(h, beta, "cf_als", c(beta = 1, h = 1),
                    call("hrfals_fit"), fmrihrf::HRF_SPMG1, "term",
                    phi, dinfo, resids)
  expect_warning(pred <- predict(fit), "Fitted values computation")
  expect_equal(pred, -resids)
})
