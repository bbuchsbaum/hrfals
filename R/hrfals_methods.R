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
#' @param recon_hrf Matrix of reconstructed HRF shapes.
#' @param gof Numeric vector of goodness-of-fit statistics per voxel.
#' @return An object of class `hrfals_fit`.
#' @export
hrfals_fit <- function(h_coeffs, beta_amps, method, lambdas, call,
                       fmrireg_hrf_basis_used, design_info, residuals,
                       recon_hrf = NULL, gof = NULL) {
  out <- list(h_coeffs = h_coeffs,
              beta_amps = beta_amps,
              method_used = method,
              lambdas = lambdas,
              call = call,
              fmrireg_hrf_basis_used = fmrireg_hrf_basis_used,
              design_info = design_info,
              residuals = residuals,
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
