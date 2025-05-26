context("cfals design helpers")

test_that("reconstruction_matrix works", {
  sf <- sampling_frame(10, TR = 1)
  phi <- reconstruction_matrix(HRF_SPMG1, sf)
  expect_equal(ncol(phi), nbasis(HRF_SPMG1))
  expect_gt(nrow(phi), 1)
})

test_that("penalty_matrix defaults to identity", {
  Rm <- penalty_matrix(HRF_SPMG1)
  expect_equal(Rm, diag(nbasis(HRF_SPMG1)))
})

test_that("convolve_timeseries_with_single_basis behaves like impulse response", {
  sf <- sampling_frame(10, TR = 1)
  raw_ts <- c(1, rep(0, 9))
  conv <- convolve_timeseries_with_single_basis(raw_ts, HRF_SPMG3, 2, sf)
  grid <- seq(0, attr(HRF_SPMG3, "span"), by = sf$TR[1])
  phi <- evaluate(HRF_SPMG3, grid)
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
  sf <- sampling_frame(20, TR = 1)
  events <- data.frame(onset = c(2, 6, 12),
                       condition = factor(c("A", "B", "A")),
                       block = 1)
  emod <- event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  des <- create_fmri_design(emod, HRF_SPMG1)
  expect_type(des, "list")
  expect_equal(length(des$X_list), 2)
  expect_equal(des$d, nbasis(HRF_SPMG1))
  expect_equal(des$k, length(des$X_list))
  expect_true(is.matrix(des$Phi))
  expect_true(is.numeric(des$h_ref_shape_norm))
})


test_that("prepare_cfals_inputs_from_fmrireg_term returns expected structure", {
  sf <- sampling_frame(20, TR = 1)
  events <- data.frame(onset = c(2, 6, 12),
                       condition = factor(c("A", "B", "A")),
                       block = 1)
  emod <- event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  Y <- matrix(rnorm(20 * 2), 20, 2)
  Z <- matrix(rnorm(20), ncol = 1)
  res <- prepare_cfals_inputs_from_fmrireg_term(Y, emod, "hrf(condition)",
                                                HRF_SPMG1, Z,
                                                hrf_shape_duration_sec = 16,
                                                hrf_shape_sample_res_sec = 1)
  expect_type(res, "list")
  expect_equal(res$d_basis_dim, nbasis(HRF_SPMG1))
  expect_equal(res$k_conditions, 2)
  expect_equal(res$n_timepoints, 20)
  expect_equal(res$v_voxels, 2)
  expect_length(res$X_list_proj, 2)
  expect_true(is.matrix(res$Phi_recon_matrix))
  expect_true(is.numeric(res$h_ref_shape_canonical_p_dim))
})

test_that("prepare_cfals_inputs_from_fmrireg_term zeroes NA rows", {
  sf <- sampling_frame(10, TR = 1)
  events <- data.frame(onset = c(1, 5),
                       condition = factor(c("A", "A")),
                       block = 1)
  emod <- event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  Y <- matrix(rnorm(10 * 1), 10, 1)
  Y[3, ] <- NA
  res <- prepare_cfals_inputs_from_fmrireg_term(Y, emod, "hrf(condition)",
                                                HRF_SPMG1)
  expect_equal(res$bad_row_idx, 3)
  expect_true(all(res$Y_proj[3, ] == 0))
  for (Xc in res$X_list_proj) {
    expect_true(all(Xc[3, ] == 0))
  }
})

