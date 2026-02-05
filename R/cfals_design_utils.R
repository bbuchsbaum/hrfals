#' CFALS Design Utilities
#'
#' Helper functions for interfacing the CF-ALS engine with the
#' fmridesign HRF basis system.
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
  vals <- fmrihrf::evaluate(hrf, grid)
  if (is.vector(vals)) matrix(vals, ncol = 1L) else as.matrix(vals)
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
#' `fmrihrf::HRF_SPMG1` and sampled on the same grid as `Phi`.
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
  # Get the number of basis functions based on the HRF object type
  d <- if (inherits(hrf_basis, "HRF")) {
    # For HRF objects from fmrihrf package
    fmrihrf::nbasis(hrf_basis)
  } else {
    stop("Unsupported HRF basis type")
  }
  
  # Use the existing design matrix from the event model
  # The issue is that the event model was created with a default HRF basis
  # We need to reconstruct with the desired basis
  
  # Extract event information and rebuild the design with the correct basis
  sample_times <- fmridesign::samples(sframe, global = TRUE)
  if (!is.numeric(sample_times) || length(sample_times) < 1 || any(!is.finite(sample_times))) {
    stop("Invalid sampling_frame: sample times must be a non-empty finite numeric vector")
  }
  if (is.unsorted(sample_times)) {
    stop("Invalid sampling_frame: sample times must be sorted in increasing order")
  }
  onset_min <- min(sample_times)
  onset_max <- max(sample_times)
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
        cond_col <- term$event_table[[var_name]]
        conditions <- if (is.factor(cond_col)) {
          levels(droplevels(cond_col))
        } else {
          unique(as.character(cond_col))
        }
        
        # For each condition, create a design matrix with d columns (one per basis function)
        for (cond in conditions) {
          # Get events for this condition
          cond_mask <- term$event_table[[var_name]] == cond
          cond_onsets <- term$onsets[cond_mask]
          
          # Create design matrix for this condition with d columns
          X_cond <- matrix(0, length(sample_times), d)
          
          # Create a timeseries with impulses at event onsets (shared across basis functions)
          ts <- numeric(length(sample_times))
          dropped_outside <- 0L
          dropped_nonfinite <- 0L
          for (onset in cond_onsets) {
            if (!is.finite(onset)) {
              dropped_nonfinite <- dropped_nonfinite + 1L
              next
            }
            if (onset < onset_min || onset > onset_max) {
              dropped_outside <- dropped_outside + 1L
              next
            }
            onset_idx <- which.min(abs(sample_times - onset))
            ts[onset_idx] <- ts[onset_idx] + 1
          }
          if (dropped_outside > 0L || dropped_nonfinite > 0L) {
            warning(sprintf(
              "Dropped %d out-of-range and %d non-finite onsets for term '%s', condition '%s'.",
              dropped_outside, dropped_nonfinite, var_name, cond
            ), call. = FALSE)
          }

          for (j in seq_len(d)) {
            X_cond[, j] <- convolve_timeseries_with_single_basis(ts, hrf_basis, j, sframe)
          }
          
          X_list[[paste0(var_name, cond)]] <- X_cond
        }
      }
    }
  }

  Phi <- reconstruction_matrix(hrf_basis, sframe)
  time_points <- seq(0, attr(hrf_basis, "span"), by = sframe$TR[1])
  h_ref <- fmrihrf::evaluate(fmrihrf::HRF_SPMG1, time_points)
  h_ref <- drop(h_ref)
  h_ref <- h_ref / max(abs(h_ref))
  if (length(h_ref) != nrow(Phi)) {
    stop("Canonical HRF length does not match Phi reconstruction grid")
  }

  list(X_list = X_list,
       d = d,
       k = length(X_list),
       Phi = Phi,
       h_ref_shape_norm = h_ref)
}
