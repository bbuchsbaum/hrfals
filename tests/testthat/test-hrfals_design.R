context("hrfals_design wrapper")

library(fmrireg)

# hrfals_design should match create_cfals_design when given the same inputs

test_that("hrfals_design replicates create_cfals_design", {
  sf <- sampling_frame(10, TR = 1)
  events <- data.frame(onset = c(2, 5),
                       condition = factor(c("A", "B")),
                       block = 1)
  Y <- matrix(rnorm(20), 10, 2)

  emod <- event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)

  direct <- create_cfals_design(Y, emod, HRF_SPMG2)
  wrap <- hrfals_design(events, TR = 1, basis = HRF_SPMG2,
                        fmri_data_obj = Y)

  expect_equal(wrap$Y_proj, direct$Y_proj)
  expect_equal(length(wrap$X_list_proj), length(direct$X_list_proj))
  expect_equal(wrap$Phi_recon_matrix, direct$Phi_recon_matrix)
})

# invalid data input should raise an error

test_that("hrfals_design validates input", {
  events <- data.frame(onset = c(1, 2), condition = factor(c("A", "B")))
  expect_error(hrfals_design(events, TR = 1, basis = HRF_SPMG2,
                             fmri_data_obj = "bad"),
               "'fmri_data_obj' must be")
})
