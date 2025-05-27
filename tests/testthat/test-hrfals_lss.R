context("hrfals_lss wrapper")

library(fmrireg)

test_that("hrfals_lss runs in both modes", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  fit <- fmrireg_cfals(dat$Y, dat$event_model, HRF_SPMG3,
                       method = "ls_svd_only")
  B_shared <- hrfals_lss(fit, dat$event_model, fmri_data_obj = dat$Y,
                          mode = "shared")
  expect_equal(nrow(B_shared), length(dat$X_list))
  expect_equal(ncol(B_shared), ncol(dat$Y))

  B_auto <- hrfals_lss(fit, dat$event_model, fmri_data_obj = dat$Y,
                        mode = "auto")
  expect_equal(dim(B_auto), c(length(dat$X_list), ncol(dat$Y)))
})
