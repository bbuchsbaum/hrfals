context("CF-ALS benchmark performance")

library(fmrireg)

# load helper functions for simulation if available
if (exists("source_test_helpers")) testthat::source_test_helpers()

benchmark_cfals <- function() {
  n_cond   <- 3
  n_trials <- 20
  TR       <- 2

  design <- simulate_simple_dataset(ncond = n_cond,
                                    nreps  = n_trials,
                                    TR     = TR,
                                    snr    = 1)

  v        <- 5000
  R_true   <- 2
  set.seed(1)
  S_true   <- matrix(rnorm(v * R_true), v, R_true)
  S_true   <- sweep(S_true, 2, sqrt(colSums(S_true^2)), "/")

  B_true <- matrix(c(1,  0.5,
                     0.6,0.4,
                     0.2,0.8), n_cond, R_true, byrow = TRUE)
  Z_true <- matrix(rnorm(n_trials * n_cond * R_true, 0, 0.2),
                   n_trials * n_cond, R_true)

  Fc <- design$clean$mat[,-1] * 0
  Ft <- Fc * 0
  on   <- design$onsets
  cond_id <- design$conditions
  for (i in seq_along(on)) {
    t_idx <- round(on[i] / TR) + 1
    Fc[t_idx, cond_id[i]] <- 1
    Ft[t_idx, i] <- 1
  }

  G_neural <- Fc %*% B_true + Ft %*% Z_true

  hrf_A <- HRF_SPMG1
  hrf_B <- HRF_BSPLINE4
  theta_vec <- rep(0, v)
  theta_vec[sample(v, 0.3 * v)] <- 1
  H_fun <- function(h_idx) if (h_idx == 0) hrf_A else hrf_B

  timegrid <- design$clean$mat[,1]
  Hmat <- sapply(theta_vec, function(idx)
    evaluate(regressor(onsets = 0, duration = 0, hrf = H_fun(idx)), timegrid))

  Y_clean <- matrix(0, length(timegrid), v)
  for (r in seq_len(R_true)) {
    Y_clean <- Y_clean +
      (convolve_design(G_neural[, r], timegrid, H_fun)) %*%
      t(S_true[, r, drop = FALSE])
  }

  eps <- replicate(v,
           simulate_noise_vector(nrow(Y_clean), TR = TR,
                                   ar = 0.4,
                                   sd = sd(Y_clean) * sqrt(1/0.35 - 1)))
  Y_noisy <- Y_clean + eps

  truth <- list(HRF_shape = theta_vec,
                B = B_true,
                Z = Z_true,
                S = S_true,
                Y_clean = Y_clean)

  basis_cfals <- HRF_FIR(12, TR)
  model_obj   <- event_model(onsets = on, durations = 0,
                             amplitudes = 1, condition = cond_id)

  runtime <- system.time({
    cf_fit <- estimate_hrf_cfals(
      matrix_dataset(Y_noisy, TR = TR),
      model_obj,
      target_event_term_name = "all",
      hrf_basis_for_cfals    = basis_cfals,
      method    = "ls_svd_1als",
      lambda_b  = 20,
      lambda_h  = 1,
      fullXtX   = TRUE)
  })

  glm_fit <- fmrireg::glm_ols(matrix_dataset(Y_noisy, TR = TR),
                              model_obj, basis_cfals)

  lss_fit <- fmrireg::glm_lss(matrix_dataset(Y_noisy, TR = TR),
                              model_obj, basis_cfals)

  list(cf = cf_fit, glm = glm_fit, lss = lss_fit,
       Hmat = Hmat, truth = truth, runtime = runtime)
}

calculate_metrics <- function(result) {
  cf_fit <- result$cf
  truth  <- result$truth
  Hmat   <- result$Hmat

  h_rmse <- sqrt(colMeans((cf_fit$recon_hrf - Hmat)^2))
  beta_true <- truth$S %*% truth$B
  beta_r <- sapply(1:nrow(beta_true), function(c)
    cor(beta_true[, c], cf_fit$beta[c, ]))

  list(h_rmse = h_rmse, beta_r = beta_r,
       runtime = result$runtime["elapsed"])
}

test_that("CF-ALS benchmark meets expectations", {
  res <- benchmark_cfals()
  metrics <- calculate_metrics(res)

  expect_lt(median(metrics$h_rmse), 0.07 * max(res$Hmat))
  expect_true(all(metrics$beta_r > 0.9))
  expect_lt(metrics$runtime, 30)
})
