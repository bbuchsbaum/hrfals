simulate_simple_data <- function(ncond = 2,
                                 nreps = 12,
                                 TR = 1,
                                 snr = 1,
                                 hrf = fmrihrf::HRF_SPMG2,
                                 seed = 123) {
  # Lightweight simulator used in several unit tests. Generates a simple
  # multi-condition design, synthesises an HRF-shaped response per onset,
  # and adds Gaussian noise at the desired SNR.
  set.seed(seed)

  total_time <- nreps * ncond * 8
  sframe <- fmridesign::sampling_frame(blocklens = total_time, TR = TR)

  onsets <- sort(runif(nreps * ncond, 0, total_time - 10))
  conditions <- rep(paste0("Cond", seq_len(ncond)), each = nreps)

  n_timepoints <- length(fmridesign::samples(sframe, global = TRUE))
  nvox <- 2
  Y_clean <- matrix(0, n_timepoints, nvox)

  for (i in seq_along(onsets)) {
    onset_idx <- round(onsets[i] / TR) + 1
    if (onset_idx <= n_timepoints - 1) {
      hrf_vals <- drop(fmrihrf::evaluate(hrf, seq(0, 20, by = TR)))
      hrf_vals <- hrf_vals / max(abs(hrf_vals))
      len <- min(length(hrf_vals), n_timepoints - onset_idx + 1)
      amplitude <- 0.5 + as.integer(factor(conditions[i])) * 0.3
      Y_clean[onset_idx:(onset_idx + len - 1), ] <- Y_clean[onset_idx:(onset_idx + len - 1), ] +
        matrix(rep(hrf_vals[seq_len(len)], nvox), len, nvox) * amplitude
    }
  }

  signal_power <- mean(Y_clean^2)
  noise_sd <- if (signal_power > 0) sqrt(signal_power / snr) else 0.1
  Y_noisy <- Y_clean + matrix(rnorm(length(Y_clean), sd = noise_sd), n_timepoints, nvox)

  list(
    noisy = Y_noisy,
    clean = Y_clean,
    onsets = onsets,
    conditions = conditions,
    sampling_frame = sframe
  )
}

simulate_noise_vector <- function(n, TR = 1, ar = 0.3, sd = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ar_coeff <- if (length(ar) == 1 && ar == 0) NULL else list(ar = ar)
  series <- stats::arima.sim(model = ar_coeff, n = n, sd = sd)
  drift <- seq_len(n) * (sd * 0.01 * TR)
  as.numeric(series + drift)
}

simulate_simple_dataset <- function(ncond = 2,
                                    nreps = 12,
                                    TR = 1,
                                    snr = 1,
                                    seed = 123) {
  # Minimal stand-in for the old fmrireg helper used by benchmark tests.
  set.seed(seed)

  total_time <- ncond * nreps * 8
  timegrid <- seq(0, total_time, by = TR)
  n_timepoints <- length(timegrid)

  onsets <- sort(runif(ncond * nreps, 0, total_time - 10))
  conditions <- rep(paste0("Cond", seq_len(ncond)), each = nreps)

  clean_mat <- matrix(0, n_timepoints, ncond + 1)
  clean_mat[, 1] <- timegrid
  for (i in seq_along(onsets)) {
    idx <- round(onsets[i] / TR) + 1
    idx <- min(idx, n_timepoints)
    cond_idx <- as.integer(factor(conditions[i])) + 1
    clean_mat[idx, cond_idx] <- clean_mat[idx, cond_idx] + 1
  }

  noise_sd <- 1 / max(snr, .Machine$double.eps)
  noisy_mat <- clean_mat
  noisy_mat[, -1] <- clean_mat[, -1] + matrix(rnorm((ncond) * n_timepoints, sd = noise_sd),
                                              n_timepoints, ncond)

  list(
    clean = list(mat = clean_mat),
    noisy = list(mat = noisy_mat),
    onsets = onsets,
    conditions = conditions,
    sampling_frame = fmridesign::sampling_frame(blocklens = total_time, TR = TR)
  )
}
