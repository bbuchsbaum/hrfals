#' Convolve timeseries with a single HRF basis function
#'
#' Helper that convolves a raw onset vector with the specified basis
#' function from an HRF object.
#'
#' @param raw_timeseries Numeric vector of length n.
#' @param hrf An object of class `HRF`.
#' @param basis_function_index Integer index of the basis function.
#' @param sampling_frame A `sampling_frame` describing the sampling grid.
#' @return Numeric vector of length n containing the convolved result.
#' @keywords internal
convolve_timeseries_with_single_basis <- function(raw_timeseries,
                                                  hrf,
                                                  basis_function_index,
                                                  sampling_frame) {
  stopifnot(inherits(hrf, "HRF"))
  nb <- fmrihrf::nbasis(hrf)
  if (basis_function_index < 1 || basis_function_index > nb) {
    stop("basis_function_index out of range")
  }
  grid <- seq(0, attr(hrf, "span"), by = sampling_frame$TR[1])
  vals <- fmrihrf::evaluate(hrf, grid)
  if (is.vector(vals)) {
    vals <- matrix(vals, ncol = 1L)
  }
  kern <- vals[, basis_function_index]
  conv_full <- stats::convolve(raw_timeseries, rev(kern), type = "open")
  conv_full[seq_along(raw_timeseries)]
}

#' Create CFALS Design Matrices from fmrireg Objects
#'
#' This function leverages fmrireg's built-in design matrix creation and
#' HRF evaluation functionality to prepare inputs for CF-ALS estimation.
#' It uses the existing `create_fmri_design` function and adds CFALS-specific
#' processing.
#'
#' @param fmri_data_obj An `fmri_dataset` or numeric matrix of BOLD data.
#' @param event_model An `event_model` object from fmrireg.
#' @param hrf_basis An HRF basis object with `nbasis > 1`.
#' @param confound_obj Optional confound matrix.
#' @param baseline_model Optional baseline model whose design matrix will
#'   be projected along with `confound_obj`.
#' @param hrf_shape_duration_sec Duration for the HRF reconstruction grid.
#' @param hrf_shape_sample_res_sec Sampling resolution for the HRF grid.
#' @param design_control List with options controlling preprocessing of
#'   the design matrices. The `standardize_predictors` element indicates
#'   whether each predictor should be z-scored prior to estimation.
#'   When `cache_design_blocks = TRUE` the lagged predictor matrices are
#'   cached in memory if the estimated footprint is below ~32GB; otherwise
#'   they are generated on the fly with a message to the user.
#' @return List with projected design matrices, reconstruction info and
#'   metadata for CFALS engines. Rows of `fmri_data_obj` containing `NA`
#'   in any voxel are zeroed out along with the corresponding rows in the
#'   design matrices and `confound_obj` (if provided). The indices of these
#'   rows are returned as `bad_row_idx`. When predictor standardization is
#'   enabled the vectors `predictor_means` and `predictor_sds` contain the
#'   scaling factors applied to each design block.
#' @export
create_cfals_design <- function(fmri_data_obj,
                               event_model,
                               hrf_basis,
                               confound_obj = NULL,
                               baseline_model = NULL,
                               hrf_shape_duration_sec = attr(hrf_basis, "span"),
                               hrf_shape_sample_res_sec = event_model$sampling_frame$TR[1],
                               design_control = list(standardize_predictors = TRUE,
                                                     cache_design_blocks = TRUE)) {
  
  # Extract BOLD data matrix
  if (inherits(fmri_data_obj, "fmri_dataset")) {
    Y_raw <- fmrireg::get_data_matrix(fmri_data_obj)
  } else if (is.matrix(fmri_data_obj)) {
    Y_raw <- fmri_data_obj
  } else {
    stop("'fmri_data_obj' must be an 'fmri_dataset' or matrix")
  }

  n_timepoints <- nrow(Y_raw)
  v_voxels <- ncol(Y_raw)

  # optional baseline design matrix
  baseline_mat <- NULL
  if (!is.null(baseline_model)) {
    baseline_mat <- fmrireg::design_matrix(baseline_model)
    if (nrow(baseline_mat) != n_timepoints) {
      stop("baseline_model design matrix has incorrect number of rows")
    }
  }

  # Handle missing data
  bad_row_idx <- which(apply(Y_raw, 1, function(r) any(is.na(r))))
  if (length(bad_row_idx) > 0) {
    Y_raw[bad_row_idx, ] <- 0
    if (!is.null(confound_obj)) {
      confound_obj[bad_row_idx, ] <- 0
    }
    if (!is.null(baseline_mat)) {
      baseline_mat[bad_row_idx, ] <- 0
    }
  }

  # Get basis dimensions
  d_basis_dim <- fmrihrf::nbasis(hrf_basis)
  if (d_basis_dim <= 1) {
    stop("CF-ALS requires an hrf_basis with nbasis > 1")
  }

  # Use the existing create_fmri_design function
  design_info <- create_fmri_design(event_model, hrf_basis)
  
  X_list_raw <- design_info$X_list
  k_conditions <- design_info$k
  Phi_recon_matrix <- design_info$Phi
  
  # Handle missing data in design matrices
  if (length(bad_row_idx) > 0) {
    X_list_raw <- lapply(X_list_raw, function(X) {
      X[bad_row_idx, ] <- 0
      X
    })
  }
  
  # Get condition names
  cond_names <- names(X_list_raw)
  
  if (k_conditions == 0) {
    stop("No estimable conditions found in event model")
  }

  # Create canonical reference HRF using the same grid as the reconstruction matrix
  # This ensures consistency between h_ref_shape_canonical and Phi_recon_matrix
  h_ref_shape_canonical <- design_info$h_ref_shape_norm

  # combine baseline model with confounds if supplied
  if (!is.null(baseline_mat)) {
    if (is.null(confound_obj)) {
      confound_obj <- baseline_mat
    } else {
      confound_obj <- cbind(confound_obj, baseline_mat)
    }
  }

  # Project out confounds using the existing function
  proj <- project_confounds(Y_raw, X_list_raw, confound_obj)
  Y_proj <- proj$Y
  X_list_proj <- proj$X_list
  
  # Zero out bad rows again after projection (QR projection can introduce non-zero values)
  if (length(bad_row_idx) > 0) {
    Y_proj[bad_row_idx, ] <- 0
    X_list_proj <- lapply(X_list_proj, function(X) {
      X[bad_row_idx, ] <- 0
      X
    })
  }

  # Determine whether to cache the projected design blocks
  size_est <- as.numeric(k_conditions) * n_timepoints * d_basis_dim * 8
  cache_requested <- isTRUE(design_control$cache_design_blocks)
  mem_threshold <- 32 * 1024^3  # ~32GB
  cache_design_blocks_used <- cache_requested && size_est <= mem_threshold
  if (cache_requested && !cache_design_blocks_used) {
    message("Estimated memory for lagged design blocks (", size_est,
            " bytes) exceeds threshold; disabling cache_design_blocks")
  }

  predictor_means <- rep(0, length(X_list_proj))
  predictor_sds <- rep(1, length(X_list_proj))
  if (isTRUE(design_control$standardize_predictors)) {
    for (i in seq_along(X_list_proj)) {
      mu_i <- mean(X_list_proj[[i]], na.rm = TRUE)
      sd_i <- stats::sd(as.vector(X_list_proj[[i]]), na.rm = TRUE)
      if (is.na(sd_i) || sd_i < .Machine$double.eps) sd_i <- 1
      X_list_proj[[i]] <- (X_list_proj[[i]] - mu_i) / sd_i
      predictor_means[i] <- mu_i
      predictor_sds[i] <- sd_i
    }
    
    # Zero out bad rows again after standardization
    if (length(bad_row_idx) > 0) {
      X_list_proj <- lapply(X_list_proj, function(X) {
        X[bad_row_idx, ] <- 0
        X
      })
    }
  }

  # Return comprehensive design information
  list(
    Y_proj = Y_proj,
    X_list_proj = X_list_proj,
    d_basis_dim = d_basis_dim,
    k_conditions = k_conditions,
    Phi_recon_matrix = Phi_recon_matrix,
    h_ref_shape_canonical = h_ref_shape_canonical,
    h_ref_shape_norm = design_info$h_ref_shape_norm,
    n_timepoints = n_timepoints,
    v_voxels = v_voxels,
    bad_row_idx = bad_row_idx,
    condition_names = cond_names,
    hrf_basis = hrf_basis,
    event_model = event_model,
    sampling_frame = event_model$sampling_frame,
    X_list = X_list_raw,
    d = d_basis_dim,
    k = k_conditions,
    Phi = Phi_recon_matrix,
    predictor_means = predictor_means,
    predictor_sds = predictor_sds,
    cache_design_blocks = cache_design_blocks_used
  )
}

#' Create CFALS Design from fmrireg Model
#'
#' @description
#' This helper previously attempted to build CFALS design matrices directly
#' from an `fmri_model` object. The approach proved brittle and has been
#' removed.  Users should instead call [create_cfals_design()] with an
#' `event_model` and HRF basis.
#'
#' @keywords internal
create_cfals_design_from_model <- function(...) {
  .Defunct("create_cfals_design",
          msg = "'create_cfals_design_from_model' has been removed.\n",
          package = "hrfals")
}

#' Legacy function name for backward compatibility
#'
#' @param ... Arguments passed to create_cfals_design
#' @keywords internal
prepare_cfals_inputs_from_fmrireg_term <- function(...) {
  .Deprecated("create_cfals_design", 
              msg = "prepare_cfals_inputs_from_fmrireg_term is deprecated. Use create_cfals_design instead.")
  create_cfals_design(...)
}

