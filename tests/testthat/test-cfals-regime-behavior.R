library(testthat)
library(fmridesign)

context("CF-ALS regime behavior")

simulate_high_k_trialwise <- function(seed) {
  set.seed(seed)

  TR <- 1
  n_time <- 220
  sf <- sampling_frame(blocklens = n_time, TR = TR)

  n_trials <- 35
  base_onsets <- seq(6, by = 6, length.out = n_trials)
  onsets <- pmin(
    base_onsets + sample(c(0, 1, 2, 3), n_trials, replace = TRUE),
    n_time - 15
  )
  onsets <- sort(unique(onsets))
  n_trials <- length(onsets)

  events <- data.frame(
    onset = onsets,
    trial_id = factor(sprintf("t%02d", seq_len(n_trials))),
    stim = factor("all"),
    block = 1
  )

  emod_trial <- event_model(
    onset ~ hrf(trial_id),
    data = events,
    block = ~ block,
    sampling_frame = sf
  )
  emod_pooled <- event_model(
    onset ~ hrf(stim),
    data = events,
    block = ~ block,
    sampling_frame = sf
  )

  basis <- fmrihrf::HRF_FIR
  design_trial <- hrfals::create_fmri_design(emod_trial, basis)
  d <- ncol(design_trial$X_list[[1]])
  k <- length(design_trial$X_list)
  v <- 20

  h_true <- c(1.3, 0.8, 0.2, rep(0, d - 3))
  beta_true <- matrix(rnorm(k * v, sd = 0.7), k, v)
  dimnames(beta_true) <- list(names(design_trial$X_list), paste0("vox", seq_len(v)))

  Y <- matrix(0, n_time, v)
  colnames(Y) <- colnames(beta_true)
  for (j in seq_len(k)) {
    signal_j <- drop(design_trial$X_list[[j]] %*% h_true)
    Y <- Y +
      matrix(rep(signal_j, v), n_time, v) *
        matrix(rep(beta_true[j, ], each = n_time), n_time, v)
  }
  Y <- Y + matrix(rnorm(length(Y), sd = 1.2), n_time, v)

  # Shared-HRF first, then LSS for trial betas.
  fit_hrf <- hrfals::hrfals(
    fmri_data_obj = Y,
    event_model = emod_pooled,
    hrf_basis = basis,
    method = "ls_svd_only",
    lam_beta = 2,
    lam_h = 2,
    max_alt = 1
  )

  # Direct trial-wise CF-ALS with fullXtX selected automatically.
  fit_direct <- hrfals::hrfals(
    fmri_data_obj = Y,
    event_model = emod_trial,
    hrf_basis = basis,
    method = "cf_als",
    lam_beta = 2,
    lam_h = 2,
    max_alt = 1,
    fullXtX = "auto"
  )

  lss <- hrfals::hrfals_lss(
    cf_fit = fit_hrf,
    events = events,
    fmri_data_obj = Y,
    formula = onset ~ hrf(trial_id),
    TR = TR,
    mode = "shared"
  )

  cor_by_voxel <- function(est, truth) {
    stopifnot(nrow(est) == nrow(truth), ncol(est) == ncol(truth))
    sapply(seq_len(ncol(truth)), function(i) cor(est[, i], truth[, i]))
  }

  direct_est <- fit_direct$beta_amps
  lss_est <- lss$betas[rownames(beta_true), colnames(beta_true), drop = FALSE]

  direct_cor <- cor_by_voxel(direct_est, beta_true)
  lss_cor <- cor_by_voxel(lss_est, beta_true)

  c(
    k = k,
    direct_median_cor = median(direct_cor, na.rm = TRUE),
    lss_median_cor = median(lss_cor, na.rm = TRUE),
    fullXtX_auto = as.numeric(fit_direct$design_info$fullXtX)
  )
}

simulate_low_k_overlap <- function(seed) {
  set.seed(seed)

  TR <- 1
  n_time <- 180
  sf <- sampling_frame(blocklens = n_time, TR = TR)

  n_events <- 14
  on1 <- seq(10, by = 10, length.out = n_events)
  on2 <- on1 + 1
  events <- data.frame(
    onset = c(on1, on2),
    condition = factor(c(rep("A", n_events), rep("B", n_events))),
    block = 1
  )
  events <- events[order(events$onset), ]

  emod <- event_model(
    onset ~ hrf(condition),
    data = events,
    block = ~ block,
    sampling_frame = sf
  )

  basis <- fmrihrf::HRF_FIR
  design <- hrfals::create_fmri_design(emod, basis)
  d <- design$d
  v <- 12

  h_true <- c(1.2, 0.8, 0.2, rep(0, d - 3))
  beta_true <- rbind(
    rnorm(v, mean = 1.2, sd = 0.15),
    rnorm(v, mean = 0.8, sd = 0.15)
  )

  Y <- matrix(0, n_time, v)
  for (cnd in 1:2) {
    signal_c <- drop(design$X_list[[cnd]] %*% h_true)
    Y <- Y +
      matrix(rep(signal_c, v), n_time, v) *
        matrix(rep(beta_true[cnd, ], each = n_time), n_time, v)
  }
  Y <- Y + matrix(rnorm(length(Y), sd = 0.15), n_time, v)

  fit_full <- hrfals::hrfals(
    fmri_data_obj = Y,
    event_model = emod,
    hrf_basis = basis,
    method = "cf_als",
    lam_beta = 0.5,
    lam_h = 0.5,
    max_alt = 1,
    fullXtX = TRUE
  )
  fit_diag <- hrfals::hrfals(
    fmri_data_obj = Y,
    event_model = emod,
    hrf_basis = basis,
    method = "cf_als",
    lam_beta = 0.5,
    lam_h = 0.5,
    max_alt = 1,
    fullXtX = FALSE
  )

  beta_abs_err <- function(est) median(abs(est - beta_true))

  c(
    full_err = beta_abs_err(fit_full$beta_amps),
    diag_err = beta_abs_err(fit_diag$beta_amps),
    full_flag = as.numeric(fit_full$design_info$fullXtX),
    diag_flag = as.numeric(fit_diag$design_info$fullXtX)
  )
}

test_that("High-k trial-wise regime favors LSS on average and auto disables fullXtX", {
  seeds <- 1:10
  res <- t(vapply(seeds, simulate_high_k_trialwise, numeric(4)))

  expect_true(all(is.finite(res)))
  expect_true(all(res[, "k"] > 12))
  expect_true(all(res[, "fullXtX_auto"] == 0))

  improvement <- res[, "lss_median_cor"] - res[, "direct_median_cor"]
  expect_gt(mean(improvement), 0.03)
  expect_gte(mean(improvement > 0), 0.60)
})

test_that("Low-k overlapping regime benefits from full cross-condition terms", {
  seeds <- 1:20
  res <- t(vapply(seeds, simulate_low_k_overlap, numeric(4)))

  expect_true(all(is.finite(res)))
  expect_true(all(res[, "full_flag"] == 1))
  expect_true(all(res[, "diag_flag"] == 0))

  improvement <- res[, "diag_err"] - res[, "full_err"]
  expect_gt(mean(improvement), 0.50)
  expect_gte(mean(improvement > 0), 0.95)
})
