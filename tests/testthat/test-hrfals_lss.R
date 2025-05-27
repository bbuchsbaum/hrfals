context("hrfals_lss wrapper")

library(fmrireg)

test_that("hrfals_lss runs in both modes", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  fit <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", HRF_SPMG3,
                           method = "ls_svd_only")
  B_shared <- hrfals_lss(fit, dat$event_model, fmri_data_obj = dat$Y,
                          mode = "shared")
  expect_s3_class(B_shared, "fastlss_fit")
  expect_equal(nrow(B_shared$betas), length(dat$X_list))
  expect_equal(ncol(B_shared$betas), ncol(dat$Y))

  B_auto <- hrfals_lss(fit, dat$event_model, fmri_data_obj = dat$Y,
                        mode = "auto")
  expect_equal(dim(B_auto$betas), c(length(dat$X_list), ncol(dat$Y)))
})

# Whitening support
test_that("hrfals_lss stores whitening matrix", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  fit <- estimate_hrf_cfals(dat$Y, dat$event_model, "hrf(condition)", HRF_SPMG3,
                           method = "ls_svd_only")
  n <- nrow(dat$Y)
  set.seed(3)
  W <- chol(crossprod(matrix(rnorm(n*n), n, n)))
  res <- hrfals_lss(fit, dat$event_model, fmri_data_obj = dat$Y,
                    mode = "shared", whitening_matrix = W)
  expect_s3_class(res, "fastlss_fit")
  expect_true(is.matrix(res$whitening_matrix))
  expect_equal(res$whitening_matrix, W)
})
