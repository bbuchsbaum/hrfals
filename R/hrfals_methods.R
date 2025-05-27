#' Construct an `hrfals_fit` object
#'
#' Simple constructor used by higher-level functions to package
#' the results returned by the various CFALS engines.
#'
#' @param h_coeffs Matrix of HRF basis coefficients (d x v).
#' @param beta_amps Matrix of condition amplitudes (k x v).
#' @param method Character string indicating the estimation method.
#' @param lambdas Numeric vector of regularisation parameters.
#' @param call The matched call to the wrapper function.
#' @param fmrireg_hrf_basis_used HRF basis object supplied to the wrapper.
#' @param design_info List with design metadata (d, k, n, v, fullXtX).
#' @param residuals Residual matrix from the projected data fit.
#' @param bad_row_idx Integer vector of time points that were zeroed due to NA
#'   values.
#' @param recon_hrf Matrix of reconstructed HRF shapes.
#' @param gof Numeric vector of goodness-of-fit statistics per voxel.
#' @return An object of class `hrfals_fit`.
#' @export
hrfals_fit <- function(h_coeffs, beta_amps, method, lambdas, call,
                       fmrireg_hrf_basis_used, target_event_term_name,
                       phi_recon_matrix, design_info, residuals,
                       bad_row_idx = integer(0),
                       recon_hrf = NULL, gof = NULL) {
  out <- list(h_coeffs = h_coeffs,
              beta_amps = beta_amps,
              method_used = method,
              lambdas = lambdas,
              call = call,
              fmrireg_hrf_basis_used = fmrireg_hrf_basis_used,
              target_event_term_name = target_event_term_name,
              phi_recon_matrix = phi_recon_matrix,
              design_info = design_info,
              residuals = residuals,
              bad_row_idx = bad_row_idx,
              reconstructed_hrfs = recon_hrf,
              gof_per_voxel = gof)
  class(out) <- c("hrfals_fit", "list")
  out
}

#' @export
print.hrfals_fit <- function(x, ...) {
  cat("\nhrfals Fit\n")
  cat("===========\n")
  info <- x$design_info
  cat(sprintf("Voxels: %d\n", info$v))
  cat(sprintf("Time points: %d\n", info$n))
  cat(sprintf("Conditions: %d\n", info$k))
  cat(sprintf("Basis functions: %d\n", info$d))
  if (!is.null(x$target_event_term_name))
    cat(sprintf("Target term: %s\n", x$target_event_term_name))
  invisible(x)
}

#' @export
summary.hrfals_fit <- function(object, ...) {
  res <- list(r2 = object$gof_per_voxel,
              design = object$design_info,
              lambdas = object$lambdas)
  class(res) <- "summary.hrfals_fit"
  res
}

#' @export
print.summary.hrfals_fit <- function(x, ...) {
  cat("\nSummary of hrfals Fit\n")
  cat("=====================\n")
  cat(sprintf("Voxels: %d\n", x$design$v))
  cat(sprintf("Time points: %d\n", x$design$n))
  cat(sprintf("Conditions: %d\n", x$design$k))
  cat(sprintf("Basis functions: %d\n", x$design$d))
  cat(sprintf("Lambda beta: %.3f\n", x$lambdas["beta"]))
  cat(sprintf("Lambda h: %.3f\n", x$lambdas["h"]))
  if (!is.null(x$r2)) {
    cat(sprintf("Mean R²: %.3f\n", mean(x$r2, na.rm = TRUE)))
    cat(sprintf("R² range: [%.3f, %.3f]\n", min(x$r2, na.rm = TRUE), max(x$r2, na.rm = TRUE)))
  }
  invisible(x)
}

#' @export
plot.hrfals_fit <- function(x, vox = 1, ...) {
  if (is.null(x$phi_recon_matrix))
    stop("phi_recon_matrix not available for plotting")
  if (vox < 1 || vox > ncol(x$h_coeffs))
    stop("'vox' out of range")
  hrf <- x$phi_recon_matrix %*% x$h_coeffs[, vox]
  plot(hrf, type = "l", xlab = "Time index", ylab = "Amplitude",
       main = paste("Reconstructed HRF - voxel", vox), ...)
  invisible(hrf)
}
