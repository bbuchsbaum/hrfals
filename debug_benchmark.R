devtools::load_all()
source('benchmark_only.R')

cat("Running benchmark...\n")
res <- benchmark_cfals()

cat("CF-ALS iterations performed:", attr(res$cf$h_coeffs, "iterations"), "\n")

# Check the reconstruction matrix details
cat("\nReconstruction matrix diagnostics:\n")
cat("Phi dimensions:", dim(res$cf$phi_recon_matrix), "\n")
cat("h_coeffs dimensions:", dim(res$cf$h_coeffs), "\n")

# Check what time grid the Phi matrix corresponds to
basis_cfals <- HRF_BSPLINE
sf <- sampling_frame(blocklens = 218, TR = 2)  # Same as in benchmark
phi_manual <- reconstruction_matrix(basis_cfals, sf)
cat("Manual Phi dimensions:", dim(phi_manual), "\n")

# Check the time grid that reconstruction_matrix uses
tr <- sf$TR[1]
span <- attr(basis_cfals, "span")
recon_timegrid <- seq(0, span, by = tr)
cat("Reconstruction time grid:", recon_timegrid, "\n")
cat("Truth time grid:", res$truth$timegrid, "\n")
cat("Grids match:", identical(recon_timegrid, res$truth$timegrid), "\n")

# Manual reconstruction to verify
manual_recon <- phi_manual %*% res$cf$h_coeffs
cat("Manual reconstruction dimensions:", dim(manual_recon), "\n")
cat("Matches CF-ALS reconstructed_hrfs:", all.equal(manual_recon, res$cf$reconstructed_hrfs), "\n")

cat("\nFirst voxel manual reconstruction:", manual_recon[, 1], "\n")
cat("First voxel CF-ALS reconstruction:", res$cf$reconstructed_hrfs[, 1], "\n")
cat("First voxel truth:", res$truth$Hmat[, 1], "\n")

cat("Calculating metrics...\n")
metrics <- calculate_metrics(res)

cat("=== RESULTS ===\n")
cat("RMSE median:", metrics$h_rmse_median, "\n")
cat("Max true:", metrics$h_rmse_max_true, "\n")
cat("Threshold (20%):", 0.2 * metrics$h_rmse_max_true, "\n")
cat("Ratio:", round(100 * metrics$h_rmse_median / metrics$h_rmse_max_true, 1), "%\n")

cat("\nFirst few RMSE values:\n")
if ("reconstructed_hrfs" %in% names(res$cf) && !is.null(res$cf$reconstructed_hrfs)) {
  h_rmse <- sqrt(colMeans((res$cf$reconstructed_hrfs - res$truth$Hmat)^2))
  cat("RMSE range:", range(h_rmse), "\n")
  cat("First 10 RMSE values:", h_rmse[1:10], "\n")
}

cat("\nFirst voxel comparison:\n")
cat("CF-ALS HRF:", res$cf$reconstructed_hrfs[, 1], "\n")
cat("Truth HRF:", res$truth$Hmat[, 1], "\n")

# Let's also check the beta correlations
cat("\nBeta correlations:", metrics$beta_r_cfals, "\n") 