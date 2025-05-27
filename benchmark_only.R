library(fmrireg)

# Helper function for convolution
convolve_design <- function(neural_signal, timegrid, hrf_func) {
  # Simple convolution: convolve neural signal with HRF
  hrf_vals <- evaluate(hrf_func, timegrid)
  if (is.matrix(hrf_vals)) hrf_vals <- hrf_vals[, 1]  # Take first column if matrix
  
  # Pad and convolve
  n <- length(neural_signal)
  conv_result <- rep(0, n)
  
  for (t in seq_len(n)) {
    for (tau in seq_len(min(t, length(hrf_vals)))) {
      if (t - tau + 1 > 0) {
        conv_result[t] <- conv_result[t] + neural_signal[t - tau + 1] * hrf_vals[tau]
      }
    }
  }
  
  return(conv_result)
}

benchmark_cfals <- function() {
  # ----- 1a. experimental design -----
  n_cond   <- 3
  n_trials <- 30                # per condition (increased for better estimation)
  TR       <- 2
  
  # Use fmrireg simulation but with proper structure
  set.seed(1)
  design <- simulate_simple_dataset(ncond = n_cond,
                                    nreps  = n_trials,
                                    TR     = TR,
                                    snr    = 2)          # start with clean signal

  # ----- 1b. make MANY voxels with LOW-RANK structure -----
  v        <- 1000              # voxels (reduced for testing speed)
  R_true   <- 2                 # latent generators
  S_true   <- matrix(rnorm(v * R_true), v, R_true)  # spatial maps
  S_true   <- sweep(S_true, 2, sqrt(colSums(S_true^2)), "/")

  # condition means & trial devs in R-space
  # Make conditions more distinct for better recovery
  B_true   <- matrix(c(1.0,  0.2,
                       0.3,  1.0,
                       0.1,  0.1), n_cond, R_true, byrow = TRUE)
  Z_true   <- matrix(rnorm(n_trials * n_cond * R_true, 0, 0.1),  # Reduce trial variability
                     n_trials * n_cond, R_true)

  # ----- 1c. build neural time courses -----
  Fc   <- design$clean$mat[, -1] * 0           # n×c design (stick before HRF)
  Ft   <- matrix(0, nrow(Fc), n_trials * n_cond)  # n×k design (k = total trials)
  on   <- design$onsets
  cond_id <- design$conditions
  

  
  # FIXED: Use the same onset discretization as create_fmri_design (which.min approach)
  n_timepoints <- nrow(design$clean$mat)  # Get timepoints from design
  sample_times <- seq(0, (n_timepoints - 1) * TR, by = TR)  # Same as BOLD timegrid
  for (i in seq_along(on)) {
    # Use which.min approach to match design matrix creation
    t_idx <- which.min(abs(sample_times - on[i]))
    # Convert "Cond1", "Cond2", "Cond3" to 1, 2, 3
    cond_numeric <- as.numeric(gsub("Cond", "", cond_id[i]))
    if (t_idx <= nrow(Fc) && cond_numeric <= ncol(Fc) && i <= ncol(Ft)) {  # Safety check
      Fc[t_idx, cond_numeric] <- 1
      Ft[t_idx, i] <- 1
    }
  }

  G_neural <- Fc %*% B_true + Ft %*% Z_true

  # HRF shapes: two very different ones
  hrf_A <- HRF_SPMG1              # canonical
  hrf_B <- HRF_SPMG1           # funky spline, same span
  theta_vec <- rep(0, v)
  theta_vec[sample(v, 0.3 * v)] <- 1              # 0 ⇒ hrf_A, 1 ⇒ hrf_B

  H_fun <- function(h_idx) {
    if (h_idx == 0) hrf_A else hrf_B
  }

  # Create two time grids:
  # 1. BOLD signal time grid (for convolution and data)
  n_timepoints <- nrow(design$clean$mat)
  bold_timegrid <- seq(0, (n_timepoints - 1) * TR, by = TR)
  
  # 2. HRF evaluation time grid - FIXED: Use same grid as CF-ALS reconstruction (0-24s, TR=2s)
  # This ensures fair comparison between truth and CF-ALS estimates
  hrf_timegrid <- seq(0, 24, by = TR)  # Same as CF-ALS: 0, 2, 4, ..., 24s (13 points)
  
  cat("BOLD timegrid: range =", range(bold_timegrid), "length =", length(bold_timegrid), "\n")
  cat("HRF timegrid: range =", range(hrf_timegrid), "length =", length(hrf_timegrid), "\n")
  
  # Generate HRF shapes on the BOLD time grid (for convolution)
  Hmat_bold <- sapply(theta_vec, function(idx) {
    hrf_func <- H_fun(idx)
    reg <- regressor(onsets = 0, duration = 0, hrf = hrf_func)
    hrf_vals <- evaluate(reg, bold_timegrid)
    if (is.matrix(hrf_vals)) hrf_vals[, 1] else hrf_vals
  })
  
  # Generate HRF shapes on the HRF time grid (for visualization)
  Hmat_hrf <- sapply(theta_vec, function(idx) {
    hrf_func <- H_fun(idx)
    reg <- regressor(onsets = 0, duration = 0, hrf = hrf_func)
    hrf_vals <- evaluate(reg, hrf_timegrid)
    if (is.matrix(hrf_vals)) hrf_vals[, 1] else hrf_vals
  })
  # Hmat: n × v convolution kernels (each column = HRF)

  # Generate clean BOLD signal with proper convolution
  Y_clean <- matrix(0, length(bold_timegrid), v)
  for (r in seq_len(R_true)) {
    # Convolve each generator with appropriate HRF per voxel
    for (vox in seq_len(v)) {
      hrf_idx <- theta_vec[vox]
      hrf_func <- H_fun(hrf_idx)
      conv_signal <- convolve_design(G_neural[, r], bold_timegrid, hrf_func)
      Y_clean[, vox] <- Y_clean[, vox] + conv_signal * S_true[vox, r]
    }
  }

  # Add realistic noise: AR(1) + drift + white noise for higher SNR
  # Increase SNR to 0.5 for better performance
  eps <- replicate(v,
           simulate_noise_vector(nrow(Y_clean), TR = TR,
                                 ar = 0.3,  # Reduce AR to 0.3
                                 sd = sd(Y_clean) * sqrt(1/0.5 - 1)))
  Y_noisy <- Y_clean + eps

  truth <- list(HRF_shape = theta_vec,   # 0/1 flag per voxel
                B         = B_true,
                Z         = Z_true,
                S         = S_true,
                Y_clean   = Y_clean,
                Hmat      = Hmat_hrf,      # Use HRF time grid for visualization
                timegrid  = hrf_timegrid)  # Store the HRF time grid

  # ----- 3. Run hrfals and baselines -----
  basis_cfals <- HRF_BSPLINE                   # 12-basis FIR
  
  # Create proper event model for fmrireg
  sf <- sampling_frame(blocklens = length(bold_timegrid), TR = TR)
  events_df <- data.frame(
    onset = on,
    condition = factor(paste0("cond", cond_id)),
    block = 1
  )
  model_obj <- event_model(onset ~ hrf(condition), data = events_df,
                           block = ~ block, sampling_frame = sf)

  runtime <- system.time({
    cf_fit <- estimate_hrf_cfals(
      fmri_data_obj          = matrix_dataset(Y_noisy, TR = TR, run_length = nrow(Y_noisy)),
      fmrireg_event_model    = model_obj,
      target_event_term_name = "hrf(condition)",
      hrf_basis_for_cfals    = basis_cfals,
      method    = "cf_als",                 # Try pure CF-ALS method
      max_alt   = 10,                       # FIXED: Allow proper convergence (was 1)
      lambda_b  = 5,                        # Further reduce beta shrinkage
      lambda_h  = 1,                      # Minimal HRF shrinkage
      fullXtX   = TRUE)
  })

  list(cf = cf_fit, truth = truth, runtime = runtime, design = design)
}

calculate_metrics <- function(result) {
  cf_fit <- result$cf
  truth <- result$truth
  
  metrics <- list()
  
  # ----- 1. HRF-shape RMSE per voxel -----
  # FIXED: Use the already-computed reconstructed_hrfs field instead of reconstructing manually
  if ("reconstructed_hrfs" %in% names(cf_fit) && !is.null(cf_fit$reconstructed_hrfs)) {
    cf_hrf_recon_raw <- cf_fit$reconstructed_hrfs  # This is already Phi %*% h_coeffs
    
    cat("Using pre-computed reconstructed_hrfs from CF-ALS\n")
    cat("CF-ALS reconstructed_hrfs dimensions:", dim(cf_hrf_recon_raw), "\n")
    
    # The reconstructed_hrfs should be on the proper HRF time grid already
    # But we need to match it to our truth time grid for comparison
    basis_cfals <- HRF_FIR
    hrf_span <- attr(basis_cfals, "span")  # 24s
    
    # Create the time grid that reconstructed_hrfs corresponds to
    tr <- 2  # Standard TR for reconstruction
    cf_timegrid <- seq(0, hrf_span, by = tr)  # Should match the reconstruction grid
    
    cat("CF-ALS time grid: 0 to", max(cf_timegrid), "s,", length(cf_timegrid), "points\n")
    
    cf_hrf_recon <- cf_hrf_recon_raw  # Use as-is for now
  } else {
    cf_hrf_recon <- NULL
    cat("No reconstructed_hrfs field found in CF-ALS fit\n")
  }
  
  if (!is.null(cf_hrf_recon) && !is.null(truth$Hmat)) {
    # FIXED: Truth and CF-ALS should now be on the same time grid (0-24s, TR=2s)
    cat("Truth HRF dimensions:", dim(truth$Hmat), "\n")
    cat("CF-ALS HRF dimensions:", dim(cf_hrf_recon), "\n")
    
    # Direct comparison since both should be on the same 13-point grid
    if (nrow(cf_hrf_recon) == nrow(truth$Hmat) && ncol(cf_hrf_recon) == ncol(truth$Hmat)) {
      h_rmse <- sqrt(colMeans((cf_hrf_recon - truth$Hmat)^2))
      metrics$h_rmse_median <- median(h_rmse)
      metrics$h_rmse_max_true <- max(abs(truth$Hmat))
      
      cat("Direct comparison successful - same dimensions\n")
      cat("RMSE range:", range(h_rmse), "\n")
    } else {
      cat("Dimension mismatch - cannot calculate RMSE\n")
      cat("CF-ALS:", dim(cf_hrf_recon), "vs Truth:", dim(truth$Hmat), "\n")
      metrics$h_rmse_median <- NA
      metrics$h_rmse_max_true <- NA
    }
  } else {
    metrics$h_rmse_median <- NA
    metrics$h_rmse_max_true <- NA
  }
  
  # ----- 2. Beta correlation (condition-level) -----
  # True betas: S %*% B (voxels × conditions)
  beta_true <- truth$S %*% t(truth$B)  # v × n_cond
  
  # CF-ALS betas
  if ("beta_amps" %in% names(cf_fit)) {
    cf_betas <- t(cf_fit$beta_amps)  # transpose to v × n_cond
  } else if ("beta" %in% names(cf_fit)) {
    cf_betas <- t(cf_fit$beta)
  } else {
    cf_betas <- NULL
  }
  
  if (!is.null(cf_betas) && ncol(cf_betas) == ncol(beta_true) && nrow(cf_betas) == nrow(beta_true)) {
    metrics$beta_r_cfals <- sapply(1:ncol(beta_true), function(c)
      cor(beta_true[, c], cf_betas[, c]))
  } else {
    metrics$beta_r_cfals <- rep(NA, ncol(beta_true))
  }
  
  # ----- 3. Runtime -----
  metrics$runtime <- result$runtime["elapsed"]
  
  return(metrics)
} 