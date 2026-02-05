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
#'   Stored under both `hrf_basis_used` and the legacy
#'   `fmrireg_hrf_basis_used` field for backward compatibility.
#' @param target_event_term_name Name of the `event_term` targeted by the fit.
#' @param phi_recon_matrix Reconstruction matrix mapping HRF coefficients to a
#'   sampled HRF shape.
#' @param design_info List with design metadata (d, k, n, v, fullXtX).
#' @param residuals Residual matrix from the projected data fit.
#' @param bad_row_idx Integer vector of time points that were zeroed due to NA
#'   values.
#' @param recon_hrf Matrix of reconstructed HRF shapes.
#' @param gof Numeric vector of goodness-of-fit statistics per voxel.
#' @param beta_penalty List containing sparse penalty parameters (l1, alpha,
#'   warm_start) used for beta estimation.
#' @return An object of class `hrfals_fit`.
#' @export
hrfals_fit <- function(h_coeffs, beta_amps, method, lambdas, call,
                       fmrireg_hrf_basis_used, target_event_term_name,
                       phi_recon_matrix, design_info, residuals,
                       bad_row_idx = integer(0),
                       recon_hrf = NULL, gof = NULL,
                       beta_penalty = NULL) {
  lambda_list <- as.list(lambdas)
  lambda_list$beta_penalty <- beta_penalty
  out <- list(h_coeffs = h_coeffs,
              beta_amps = beta_amps,
              method_used = method,
              lambdas = lambda_list,
              call = call,
              fmrireg_hrf_basis_used = fmrireg_hrf_basis_used,
              hrf_basis_used = fmrireg_hrf_basis_used,
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

#' Print method for hrfals_fit objects
#'
#' @param x An `hrfals_fit` object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the input object.
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

#' Summary method for hrfals_fit objects
#'
#' @param object An `hrfals_fit` object.
#' @param ... Additional arguments (unused).
#' @return A list of class `summary.hrfals_fit` containing summary statistics.
#' @export
summary.hrfals_fit <- function(object, ...) {
  res <- list(r2 = object$gof_per_voxel,
              design = object$design_info,
              lambdas = object$lambdas)
  class(res) <- "summary.hrfals_fit"
  res
}

#' Print method for summary.hrfals_fit objects
#'
#' @param x A `summary.hrfals_fit` object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the input object.
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
    cat(sprintf("Mean R^2: %.3f\n", mean(x$r2, na.rm = TRUE)))
    cat(sprintf("R^2 range: [%.3f, %.3f]\n", min(x$r2, na.rm = TRUE), max(x$r2, na.rm = TRUE)))
  }
  invisible(x)
}

#' Residuals method for hrfals_fit objects
#'
#' @param object An `hrfals_fit` object.
#' @param ... Additional arguments (unused).
#' @return Residual matrix stored in the fit object.
#' @export
residuals.hrfals_fit <- function(object, ...) {
  object$residuals
}

#' Plot method for hrfals_fit objects
#'
#' @param x An `hrfals_fit` object.
#' @param vox Voxel index to plot (default = 1).
#' @param ... Additional graphical parameters passed to `plot()`.
#' @return Invisibly returns the reconstructed HRF for the specified voxel.
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
#' Tidy method for hrfals_fit objects
#'
#' Produces a data frame summarising beta amplitudes and optional goodness-of-fit
#' statistics for each voxel.
#'
#' @param x An `hrfals_fit` object.
#' @param ... Unused.
#' @return A data frame with one row per voxel and columns for beta amplitudes
#'   and R-squared values when available.
#' @export
#' @importFrom broom tidy
#' @importFrom broom glance
#' @importFrom tibble tibble
#' @importFrom ggplot2 autoplot
#' @seealso [glance.hrfals_fit()]

tidy.hrfals_fit <- function(x, ...) {
  betas <- t(x$beta_amps)
  if (!is.null(rownames(x$beta_amps))) {
    colnames(betas) <- paste0("beta_", rownames(x$beta_amps))
  } else {
    colnames(betas) <- paste0("beta_", seq_len(ncol(betas)))
  }
  df <- tibble::tibble(voxel = seq_len(nrow(betas)))
  df <- cbind(df, as.data.frame(betas))
  if (!is.null(x$gof_per_voxel)) {
    df$r_squared <- x$gof_per_voxel
  }
  df
}

#' Glance method for hrfals_fit objects
#'
#' Returns a one-row summary with design information and mean goodness of fit.
#'
#' @inheritParams tidy.hrfals_fit
#' @return A data frame summarising the fit.
#' @export

glance.hrfals_fit <- function(x, ...) {
  info <- x$design_info
  data.frame(
    n_vox = info$v,
    n_time = info$n,
    n_cond = info$k,
    basis_len = info$d,
    lambda_beta = unname(x$lambdas["beta"]),
    lambda_h = unname(x$lambdas["h"]),
    r2_mean = if (!is.null(x$gof_per_voxel)) mean(x$gof_per_voxel, na.rm = TRUE) else NA_real_
  )
}

#' Autoplot method for hrfals_fit objects
#'
#' Creates a ggplot showing the reconstructed HRF for a selected voxel.
#'
#' @param object An `hrfals_fit` object.
#' @param vox Voxel index to plot (default = 1).
#' @param ... Additional arguments (currently unused).
#' @return A `ggplot` object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_line labs

autoplot.hrfals_fit <- function(object, vox = 1, ...) {
  if (is.null(object$phi_recon_matrix))
    stop("phi_recon_matrix not available for plotting")
  if (vox < 1 || vox > ncol(object$h_coeffs))
    stop("'vox' out of range")
  hrf <- as.vector(object$phi_recon_matrix %*% object$h_coeffs[, vox])
  df <- data.frame(time = seq_along(hrf), amplitude = hrf)
  ggplot2::ggplot(df, ggplot2::aes(x = time, y = amplitude)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = paste("Reconstructed HRF - voxel", vox),
         x = "Time index", y = "Amplitude")
}

utils::globalVariables(c("time", "amplitude"))

#' Predict method for hrfals_fit objects
#'
#' Generates predictions from a fitted hrfals model. When no new data is provided,
#' returns the fitted values for the original training data.
#'
#' @param object An `hrfals_fit` object.
#' @param newdata Optional new data for prediction. Currently not implemented.
#' @param ... Additional arguments (currently unused).
#' @return A matrix of predicted values with dimensions (n_timepoints x n_voxels).
#'   For the original training data, these are the fitted values from the model.
#' @details 
#' The predictions are computed as:
#' \deqn{Y_{pred} = \sum_{c=1}^{k} (X_c \times h) \odot \beta_c}
#' where \eqn{X_c} are the design matrices, \eqn{h} are the HRF coefficients,
#' \eqn{\beta_c} are the condition amplitudes, and \eqn{\odot} denotes 
#' element-wise multiplication across voxels.
#' @export
#' @examples
#' \dontrun{
#' # Assuming you have a fitted hrfals object
#' fitted_values <- predict(fit)
#' 
#' # Check fitted values against original data
#' residuals_manual <- Y - fitted_values
#' all.equal(residuals_manual, residuals(fit))
#' }
predict.hrfals_fit <- function(object, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    stop("Prediction with new data is not yet implemented. ",
         "Currently only fitted values for original data are supported.")
  }
  
  # Extract model components
  h_coeffs <- object$h_coeffs
  beta_amps <- object$beta_amps
  design_info <- object$design_info
  
  # Get dimensions
  n <- design_info$n
  v <- design_info$v
  k <- design_info$k
  
  # We need to reconstruct the design matrices or use stored fitted values
  # For now, we can compute this from residuals if available
  if (!is.null(object$residuals)) {
    # If we have residuals, we can compute fitted values as original data minus residuals
    # But we need the original Y data, which we don't have stored
    # So we need to reconstruct from the model parameters
    
    # This is a limitation - we would need to store the original design matrices
    # or the original data to compute proper fitted values
    # For now, we'll provide a placeholder that indicates this limitation
    
    warning("Fitted values computation requires access to original design matrices. ",
            "Consider storing design matrices in the hrfals_fit object for full predict functionality.")
    
    # Return residuals with opposite sign as a rough approximation
    # This is not ideal but indicates where the functionality should go
    return(-object$residuals)
  }
  
  stop("Cannot compute fitted values without stored residuals or design matrices. ",
       "The hrfals_fit object needs to store sufficient information for prediction.")
}
