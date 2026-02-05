context("cfals design helpers")

# Load required packages in the correct order
# fmridesign must be loaded for hrf() and other functions
# fmrihrf must be available for HRF basis objects
library(fmridesign)
if (!requireNamespace("fmrihrf", quietly = TRUE)) {
  skip("fmrihrf package not available")
}

test_that("reconstruction_matrix works", {
  sf <- fmridesign::sampling_frame(10, TR = 1)
  phi <- hrfals::reconstruction_matrix(fmrihrf::HRF_SPMG3, sf)
  expect_equal(ncol(phi), fmrihrf::nbasis(fmrihrf::HRF_SPMG3))
  expect_gt(nrow(phi), 1)
})

test_that("fmrihrf penalty_matrix works", {
  Rm <- fmrihrf::penalty_matrix(fmrihrf::HRF_SPMG3)
  expect_true(is.matrix(Rm))
  expect_equal(dim(Rm), c(fmrihrf::nbasis(fmrihrf::HRF_SPMG3), fmrihrf::nbasis(fmrihrf::HRF_SPMG3)))
  expect_true(isSymmetric(Rm))
})

test_that("convolve_timeseries_with_single_basis behaves like impulse response", {
  sf <- fmridesign::sampling_frame(10, TR = 1)
  raw_ts <- c(1, rep(0, 9))
  conv <- hrfals:::convolve_timeseries_with_single_basis(raw_ts, fmrihrf::HRF_SPMG3, 2, sf)
  grid <- seq(0, attr(fmrihrf::HRF_SPMG3, "span"), by = sf$TR[1])
  phi <- fmrihrf::evaluate(fmrihrf::HRF_SPMG3, grid)
  if (is.vector(phi)) phi <- matrix(phi, ncol = 1L)
  expect_equal(conv, phi[seq_along(raw_ts), 2])
})

test_that("project_confounds projects via QR", {
  X <- matrix(rnorm(20), 5, 4)
  Y <- matrix(rnorm(10), 5, 2)
  Z <- matrix(seq_len(5), ncol = 1)
  res <- project_confounds(Y, list(X), Z)
  expect_equal(dim(res$X_list[[1]]), dim(X))
  expect_equal(dim(res$Y), dim(Y))
})

test_that("create_fmri_design returns expected structure", {
  sf <- fmridesign::sampling_frame(20, TR = 1)
  events <- data.frame(onset = c(2, 6, 12),
                       condition = factor(c("A", "B", "A")),
                       block = 1)
  emod <- fmridesign::event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  des <- hrfals::create_fmri_design(emod, fmrihrf::HRF_SPMG3)
  expect_type(des, "list")
  expect_equal(length(des$X_list), 2)
  expect_equal(des$d, fmrihrf::nbasis(fmrihrf::HRF_SPMG3))
  expect_equal(des$k, length(des$X_list))
  expect_true(is.matrix(des$Phi))
  expect_true(is.numeric(des$h_ref_shape_norm))
  expect_equal(length(des$h_ref_shape_norm), nrow(des$Phi))
  expect_equal(max(abs(des$h_ref_shape_norm)), 1)
})

test_that("create_cfals_design returns expected structure", {
  sf <- fmridesign::sampling_frame(20, TR = 1)
  events <- data.frame(onset = c(2, 6, 12),
                       condition = factor(c("A", "B", "A")),
                       block = 1)
  emod <- fmridesign::event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  Y <- matrix(rnorm(20 * 2), 20, 2)
  Z <- matrix(rnorm(20), ncol = 1)
  
  res <- hrfals::create_cfals_design(Y, emod, fmrihrf::HRF_SPMG3, Z,
                            hrf_shape_duration_sec = 16,
                            hrf_shape_sample_res_sec = 1)
  
  expect_type(res, "list")
  expect_equal(res$d_basis_dim, fmrihrf::nbasis(fmrihrf::HRF_SPMG3))
  expect_equal(res$k_conditions, 2)
  expect_equal(res$n_timepoints, 20)
  expect_equal(res$v_voxels, 2)
  expect_length(res$X_list_proj, 2)
  expect_true(is.matrix(res$Phi_recon_matrix))
  expect_true(is.numeric(res$h_ref_shape_canonical))
  expect_equal(length(res$h_ref_shape_canonical), nrow(res$Phi_recon_matrix))
  expect_equal(max(abs(res$h_ref_shape_canonical)), 1)
  expect_equal(length(res$condition_names), 2)
  
  # Check compatibility fields
  expect_equal(res$d, res$d_basis_dim)
  expect_equal(res$k, res$k_conditions)
  expect_equal(res$Phi, res$Phi_recon_matrix)
  expect_length(res$X_list, 2)
  expect_true(is.numeric(res$h_ref_shape_norm))
})

test_that("create_cfals_design uses hrf_shape_duration_sec and hrf_shape_sample_res_sec", {
  sf <- fmridesign::sampling_frame(20, TR = 1)
  events <- data.frame(onset = c(2, 6, 12),
                       condition = factor(c("A", "B", "A")),
                       block = 1)
  emod <- fmridesign::event_model(onset ~ hrf(condition), data = events,
                                  block = ~ block, sampling_frame = sf)
  Y <- matrix(rnorm(20 * 2), 20, 2)

  res <- hrfals::create_cfals_design(
    Y, emod, fmrihrf::HRF_SPMG3,
    hrf_shape_duration_sec = 8,
    hrf_shape_sample_res_sec = 0.5
  )

  expect_equal(nrow(res$Phi_recon_matrix), length(seq(0, 8, by = 0.5)))
  expect_equal(length(res$h_ref_shape_canonical), nrow(res$Phi_recon_matrix))
})

test_that("create_fmri_design drops out-of-range events (with warning)", {
  sf <- fmridesign::sampling_frame(10, TR = 1)
  events_bad <- data.frame(
    onset = c(2, 20),
    condition = factor(c("A", "A")),
    block = 1
  )
  events_ok <- data.frame(
    onset = 2,
    condition = factor("A"),
    block = 1
  )

  emod_bad <- fmridesign::event_model(onset ~ hrf(condition), data = events_bad,
                                      block = ~ block, sampling_frame = sf)
  emod_ok <- fmridesign::event_model(onset ~ hrf(condition), data = events_ok,
                                     block = ~ block, sampling_frame = sf)

  des_ok <- hrfals::create_fmri_design(emod_ok, fmrihrf::HRF_SPMG3)
  des_bad <- expect_warning(
    hrfals::create_fmri_design(emod_bad, fmrihrf::HRF_SPMG3),
    "Dropped"
  )

  expect_equal(des_bad$X_list[[1]], des_ok$X_list[[1]])
})

test_that("create_fmri_design accumulates multiple events at the same sample", {
  sf <- fmridesign::sampling_frame(10, TR = 1)
  events_one <- data.frame(
    onset = 2.1,
    condition = factor("A"),
    block = 1
  )
  events_two <- data.frame(
    onset = c(2.1, 2.2),
    condition = factor(c("A", "A")),
    block = 1
  )

  emod_one <- fmridesign::event_model(onset ~ hrf(condition), data = events_one,
                                      block = ~ block, sampling_frame = sf)
  emod_two <- fmridesign::event_model(onset ~ hrf(condition), data = events_two,
                                      block = ~ block, sampling_frame = sf)

  des_one <- hrfals::create_fmri_design(emod_one, fmrihrf::HRF_SPMG3)
  des_two <- hrfals::create_fmri_design(emod_two, fmrihrf::HRF_SPMG3)

  expect_equal(des_two$X_list[[1]], 2 * des_one$X_list[[1]])
})

test_that("create_cfals_design zeroes NA rows", {
  sf <- fmridesign::sampling_frame(10, TR = 1)
  events <- data.frame(onset = c(1, 5),
                       condition = factor(c("A", "A")),
                       block = 1)
  emod <- fmridesign::event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  Y <- matrix(rnorm(10 * 1), 10, 1)
  Y[3, ] <- NA

  res <- hrfals::create_cfals_design(Y, emod, fmrihrf::HRF_SPMG3)
  
  expect_equal(res$bad_row_idx, 3)
  expect_true(all(res$Y_proj[3, ] == 0))
  for (Xc in res$X_list_proj) {
    expect_true(all(Xc[3, ] == 0))
  }
})

test_that("create_cfals_design zeroes NA rows with confounds", {
  sf <- fmridesign::sampling_frame(10, TR = 1)
  events <- data.frame(onset = c(1, 5),
                       condition = factor(c("A", "A")),
                       block = 1)
  emod <- fmridesign::event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  Y <- matrix(rnorm(10 * 1), 10, 1)
  Y[4, ] <- NA
  Z <- matrix(rnorm(10), ncol = 1)

  res <- hrfals::create_cfals_design(Y, emod, fmrihrf::HRF_SPMG3, Z)

  expect_equal(res$bad_row_idx, 4)
  expect_true(all(res$Y_proj[4, ] == 0))
  for (Xc in res$X_list_proj) {
    expect_true(all(Xc[4, ] == 0))
  }
})

test_that("create_cfals_design works without confounds", {
  sf <- fmridesign::sampling_frame(15, TR = 1)
  events <- data.frame(onset = c(2, 8),
                       condition = factor(c("A", "B")),
                       block = 1)
  emod <- fmridesign::event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  Y <- matrix(rnorm(15 * 3), 15, 3)
  
  res <- hrfals::create_cfals_design(Y, emod, fmrihrf::HRF_SPMG2)
  
  expect_type(res, "list")
  expect_equal(res$d_basis_dim, fmrihrf::nbasis(fmrihrf::HRF_SPMG2))
  expect_equal(res$k_conditions, 2)
  expect_equal(res$n_timepoints, 15)
  expect_equal(res$v_voxels, 3)
  expect_length(res$X_list_proj, 2)
  expect_length(res$condition_names, 2)
  expect_true(all(c("conditionA", "conditionB") %in% res$condition_names))
})
