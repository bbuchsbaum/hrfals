#' CFALS Design Utilities
#'
#' Helper functions for interfacing the CF-ALS engine with the
#' fmrireg HRF basis system.
#'
#' @name cfals_design_utils
NULL

#' Reconstruction matrix for an HRF basis
#'
#' Returns a matrix \eqn{\Phi} that converts basis coefficients into a
#' sampled HRF shape.
#'
#' @param hrf An object of class `HRF`.
#' @param sframe A `sampling_frame` object or numeric vector of times.
#' @return A numeric matrix with one column per basis function.
#' @export
reconstruction_matrix <- function(hrf, sframe) {
  UseMethod("reconstruction_matrix")
}

#' @export
reconstruction_matrix.HRF <- function(hrf, sframe) {
  grid <- if (inherits(sframe, "sampling_frame")) {
    seq(0, attr(hrf, "span"), by = sframe$TR[1])
  } else {
    as.numeric(sframe)
  }
  vals <- evaluate(hrf, grid)
  if (is.vector(vals)) matrix(vals, ncol = 1L) else as.matrix(vals)
}



#' Convolve a timeseries with a single HRF basis function
#'
#' Utility helper that extracts one basis function from an `HRF` object and
#' performs discrete convolution with a raw timeseries.  The result has the
#' same length as the input series and is truncated to the sampling frame.
#'
#' @param ts Numeric vector of raw onset values.
#' @param hrf_basis An object of class `HRF` providing the basis set.
#' @param basis_index Integer index of the basis function to use.
#' @param sframe A `sampling_frame` object describing the temporal grid.
#' @return Numeric vector of convolved values.
#' @keywords internal
convolve_timeseries_with_single_basis <- function(ts, hrf_basis,
                                                  basis_index = 1, sframe) {
  if (!is.numeric(ts)) {
    stop("'ts' must be numeric")
  }
  if (!inherits(hrf_basis, "HRF")) {
    stop("'hrf_basis' must be an object of class 'HRF'")
  }
  nb <- fmrireg::nbasis(hrf_basis)
  if (basis_index < 1 || basis_index > nb) {
    stop("'basis_index' out of range")
  }

  grid <- if (inherits(sframe, "sampling_frame")) {
    seq(0, attr(hrf_basis, "span"), by = sframe$TR[1])
  } else {
    as.numeric(sframe)
  }
  vals <- fmrireg::evaluate(hrf_basis, grid)
  if (is.vector(vals)) vals <- matrix(vals, ncol = 1L)
  phi_j <- vals[, basis_index]

  conv_full <- stats::convolve(ts, rev(phi_j), type = "open")
  conv_full[seq_along(ts)]
}

#' Project design and data matrices to the null space of confounds
#'
#' Projects both the data matrix `Y` and each design matrix in
#' `X_list` using QR decomposition of the confound matrix.  The
#' projection can optionally use LAPACK's QR implementation for
#' improved numerical stability.
#'
#' @param Y Numeric matrix of BOLD data (time points \eqn{\times}
#'   voxels).
#' @param X_list A list of design matrices with the same number of
#'   rows as `Y`.
#' @param confounds Optional confound matrix with matching rows.
#' @param lapack_qr Logical; passed to `qr()` as the `LAPACK`
#'   argument.
#' @return A list with projected `X_list` and `Y` matrices.
#' @export
project_confounds <- function(Y, X_list, confounds = NULL, lapack_qr = FALSE) {
  if (is.null(confounds)) {
    return(list(X_list = X_list, Y = Y))
  }
  qrZ <- qr(confounds, LAPACK = lapack_qr)
  Xp <- lapply(X_list, function(X) qr.resid(qrZ, X))
  Yp <- qr.resid(qrZ, Y)
  list(X_list = Xp, Y = Yp)
}

#' Create design matrices for CFALS estimation
#'
#' Convenience helper that constructs the list of design matrices for
#' a given `event_model` and HRF basis.  It also returns useful
#' metadata such as the number of basis functions and conditions as
#' well as a reconstruction matrix for converting HRF coefficients to
#' sampled shapes and a normalised reference HRF vector for sign
#' alignment. The reference HRF is generated using
#' `fmrireg::HRF_SPMG1` and sampled on the same grid as `Phi`.
#'
#' @param event_model An object of class `event_model`.
#' @param hrf_basis An `HRF` basis object.
#' @return A list with elements `X_list`, `d`, `k`, `Phi`, and
#'   `h_ref_shape_norm`.
#' @export
create_fmri_design <- function(event_model, hrf_basis) {
  if (!inherits(event_model, "event_model")) {
    stop("'event_model' must be an 'event_model' object")
  }
  if (!inherits(hrf_basis, "HRF")) {
    stop("'hrf_basis' must be an object of class 'HRF'")
  }

  sframe <- event_model$sampling_frame
  sample_times <- samples(sframe, global = TRUE)
  d <- nbasis(hrf_basis)

  # Extract event information from the event model
  # We need to manually create design matrices for each condition and basis function
  X_list <- list()
  
  # Get the event terms from the model
  terms <- event_model$terms
  
  # For each term that involves HRF convolution
  for (term in terms) {
    if (inherits(term, "event_term")) {
      # Get the variable name for this term
      var_name <- term$varname
      
      # Get the conditions for this term from the event_table
      if (var_name %in% names(term$event_table)) {
        conditions <- levels(term$event_table[[var_name]])
        
        # For each condition, create a design matrix with d columns (one per basis function)
        for (cond in conditions) {
          # Get events for this condition
          cond_mask <- term$event_table[[var_name]] == cond
          cond_onsets <- term$onsets[cond_mask]
          
          # Create design matrix for this condition with d columns
          X_cond <- matrix(0, length(sample_times), d)
          
          # For each basis function
          for (j in seq_len(d)) {
            # Create a timeseries with impulses at event onsets
            ts <- rep(0, length(sample_times))
            for (onset in cond_onsets) {
              onset_idx <- which.min(abs(sample_times - onset))
              if (onset_idx <= length(ts)) {
                ts[onset_idx] <- 1
              }
            }
            
            # Convolve with the j-th basis function
            X_cond[, j] <- convolve_timeseries_with_single_basis(ts, hrf_basis, j, sframe)
          }
          
          X_list[[paste0(var_name, cond)]] <- X_cond
        }
      }
    }
  }

  Phi <- reconstruction_matrix(hrf_basis, sframe)
  time_points <- seq(0, attr(hrf_basis, "span"), by = sframe$TR[1])
  h_ref <- fmrireg::evaluate(fmrireg::HRF_SPMG1, time_points)
  h_ref <- drop(h_ref)
  h_ref <- h_ref / max(abs(h_ref))
  if (length(h_ref) != nrow(Phi)) {
    stop("Canonical HRF length does not match Phi reconstruction grid")
  }

  list(X_list = X_list,
       d = nbasis(hrf_basis),
       k = length(X_list),
       Phi = Phi,
       h_ref_shape_norm = h_ref)
}
