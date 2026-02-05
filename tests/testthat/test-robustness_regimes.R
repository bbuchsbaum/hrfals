test_that("create_fmri_design warns and drops out-of-range onsets", {
  sf <- fmridesign::sampling_frame(blocklens = 20, TR = 1)
  events <- data.frame(
    onset = c(5, 25), # 25s exceeds run length but is allowed by event_model
    condition = factor(c("A", "A")),
    block = 1
  )
  emod <- fmridesign::event_model(onset ~ fmridesign::hrf(condition),
                                  data = events,
                                  block = ~ block,
                                  sampling_frame = sf)

  warns <- character()
  design <- withCallingHandlers(
    create_fmri_design(emod, fmrihrf::HRF_FIR),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("Dropped [0-9]+ out-of-range", warns)))
  expect_true(is.list(design$X_list))
  expect_gt(length(design$X_list), 0)
  expect_true(all(vapply(design$X_list, nrow, integer(1)) == length(fmridesign::samples(sf, global = TRUE))))
})

test_that("estimate_hrf_cfals returns NA gof (not NaN) for constant voxels", {
  set.seed(123)
  dat <- simulate_cfals_wrapper_data(fmrihrf::HRF_SPMG3, noise_sd = 0.05)
  Y <- dat$Y
  Y[, 1] <- 3.14 # constant voxel => SST = 0

  fit <- estimate_hrf_cfals(
    fmri_data_obj = Y,
    fmridesign_event_model = dat$event_model,
    target_event_term_name = "hrf(condition)",
    hrf_basis_for_cfals = fmrihrf::HRF_SPMG3,
    method = "cf_als",
    lambda_b = 0.1,
    lambda_h = 0.1,
    max_alt = 2,
    design_control = list(standardize_predictors = FALSE)
  )

  expect_equal(length(fit$gof), ncol(Y))
  expect_true(is.na(fit$gof[1]))
  expect_false(is.nan(fit$gof[1]))
})

test_that("cf_als_engine handles highly collinear conditions without exploding", {
  sf <- fmridesign::sampling_frame(blocklens = 40, TR = 1)
  # Two conditions with identical onsets => extreme collinearity.
  # fmridesign requires strictly increasing onsets within each block, so
  # we jitter the second condition by an epsilon that still maps to the
  # same sample index given TR = 1.
  eps <- 1e-6
  events <- data.frame(
    onset = c(5, 5 + eps, 15, 15 + eps, 25, 25 + eps, 35, 35 + eps),
    condition = factor(rep(c("A", "B"), each = 4)),
    block = 1
  )
  emod <- fmridesign::event_model(onset ~ fmridesign::hrf(condition),
                                  data = events,
                                  block = ~ block,
                                  sampling_frame = sf)
  design <- create_fmri_design(emod, fmrihrf::HRF_FIR)
  n <- length(fmridesign::samples(sf, global = TRUE))
  v <- 3

  # Small synthetic Y; regularization should keep solution bounded.
  Y <- matrix(rnorm(n * v, sd = 0.2), n, v)
  res <- hrfals:::cf_als_engine(
    X_list_proj = design$X_list,
    Y_proj = Y,
    lambda_init = 1,
    lambda_b = 1,
    lambda_h = 1,
    lambda_joint = 0.1,
    Phi_recon_matrix = design$Phi,
    h_ref_shape_canonical = rep(1, nrow(design$Phi)),
    max_alt = 3,
    fullXtX_flag = TRUE,
    precompute_xty_flag = TRUE,
    design_control = list(standardize_predictors = FALSE, cache_design_blocks = TRUE)
  )

  expect_equal(dim(res$beta), c(design$k, v))
  expect_equal(dim(res$h), c(design$d, v))
  expect_true(all(is.finite(res$beta)))
  expect_true(all(is.finite(res$h)))
  obj <- attr(res$h, "objective_trace")
  expect_true(all(is.finite(obj)))
  expect_lt(max(obj) / min(obj), 1e6) # very loose guard against numerical blow-ups
})

test_that("CG h-update warns on non-convergence but returns finite outputs", {
  dat <- simulate_cfals_wrapper_data(fmrihrf::HRF_SPMG3, noise_sd = 0.1)
  prep <- create_cfals_design(dat$Y, dat$event_model, fmrihrf::HRF_SPMG3,
                              design_control = list(standardize_predictors = FALSE))
  # line graph with one node per voxel
  v <- ncol(prep$Y_proj)
  mask <- array(1, dim = c(v, 1, 1))
  lap_obj <- build_voxel_laplacian(mask)

  warns <- character()
  res <- withCallingHandlers(
    hrfals:::cf_als_engine(
      X_list_proj = prep$X_list_proj,
      Y_proj = prep$Y_proj,
      lambda_init = 1,
      lambda_b = 0.1,
      lambda_h = 0.1,
      lambda_s = 0.1,
      laplacian_obj = lap_obj,
      h_solver = "cg",
      cg_max_iter = 1,
      cg_tol = 1e-12,
      Phi_recon_matrix = prep$Phi_recon_matrix,
      h_ref_shape_canonical = prep$h_ref_shape_canonical,
      max_alt = 1,
      design_control = list(standardize_predictors = FALSE, cache_design_blocks = TRUE)
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("CG solver did not converge", warns)))
  expect_true(all(is.finite(res$beta)))
  expect_true(all(is.finite(res$h)))
})

test_that("Sparse beta step falls back to ridge when glmnet errors", {
  dat <- simulate_cfals_wrapper_data(fmrihrf::HRF_SPMG3, noise_sd = 0.1)
  prep <- create_cfals_design(dat$Y, dat$event_model, fmrihrf::HRF_SPMG3,
                              design_control = list(standardize_predictors = FALSE))

  testthat::local_mocked_bindings(
    glmnet = function(...) stop("boom"),
    .package = "glmnet"
  )

  warns <- character()
  res <- withCallingHandlers(
    hrfals:::cf_als_engine(
      X_list_proj = prep$X_list_proj,
      Y_proj = prep$Y_proj,
      lambda_init = 1,
      lambda_b = 0.1,
      lambda_h = 0.1,
      lambda_joint = 0,
      Phi_recon_matrix = prep$Phi_recon_matrix,
      h_ref_shape_canonical = prep$h_ref_shape_canonical,
      max_alt = 1,
      beta_penalty = list(l1 = 0.01, alpha = 1, warm_start = TRUE),
      design_control = list(standardize_predictors = FALSE, cache_design_blocks = TRUE)
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("glmnet failed", warns)))
  expect_true(all(is.finite(res$beta)))
  expect_true(all(is.finite(res$h)))
})
