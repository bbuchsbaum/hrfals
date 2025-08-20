context("baseline_model handling")

library(fmrireg)

test_that("baseline_model is projected with confounds", {
  sf <- sampling_frame(10, TR = 1)
  events <- data.frame(onset = c(2, 7),
                       condition = factor(c("A", "B")),
                       block = 1)
  emod <- event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  Y <- matrix(rnorm(20), 10, 2)
  Z <- matrix(rnorm(10), ncol = 1)
  bmod <- baseline_model(sframe = sf, degree = 1)

  res1 <- create_cfals_design(Y, emod, fmrihrf::HRF_SPMG2,
                              confound_obj = Z,
                              baseline_model = bmod)
  res2 <- create_cfals_design(Y, emod, fmrihrf::HRF_SPMG2,
                              confound_obj = cbind(Z, design_matrix(bmod)))

  expect_equal(res1$Y_proj, res2$Y_proj)
  for (i in seq_along(res1$X_list_proj)) {
    expect_equal(res1$X_list_proj[[i]], res2$X_list_proj[[i]])
  }
})
