simulate_simple_data <- function(ncond = 2,
                                 nreps = 12,
                                 TR = 1,
                                 snr = 1,
                                 hrf = fmrihrf::HRF_SPMG2,
                                 seed = 123) {
  # Custom simulation function to work around fmrireg bug
  set.seed(seed)
  
  # Create sampling frame
  total_time <- nreps * ncond * 8  # 8 seconds per trial
  sframe <- fmrireg::sampling_frame(blocklens = total_time, TR = TR)
  
  # Generate random onsets
  onsets <- sort(runif(nreps * ncond, 0, total_time - 10))
  conditions <- rep(paste0("Cond", 1:ncond), each = nreps)
  
  # Simulate simple data without complex HRF modeling
  n_timepoints <- length(fmrireg::samples(sframe, global = TRUE))
  nvox <- 2  # 2 voxels
  
  # Create simple signal with some structure
  Y_clean <- matrix(0, n_timepoints, nvox)
  for (i in seq_along(onsets)) {
    onset_idx <- round(onsets[i] / TR) + 1
    if (onset_idx <= n_timepoints - 10) {
      # Add a simple response pattern
      response <- exp(-((1:10) - 3)^2 / 4)  # Simple HRF-like response
      Y_clean[onset_idx:(onset_idx + 9), ] <- Y_clean[onset_idx:(onset_idx + 9), ] + 
        matrix(rep(response, nvox), 10, nvox) * runif(1, 0.5, 1.5)
    }
  }
  
  # Add noise
  noise_sd <- sqrt(mean(Y_clean^2) / snr)
  Y_noisy <- Y_clean + matrix(rnorm(length(Y_clean), sd = noise_sd), 
                              n_timepoints, nvox)
  
  list(
    noisy = Y_noisy,
    clean = Y_clean,
    onsets = onsets,
    conditions = conditions,
    sampling_frame = sframe
  )
}
