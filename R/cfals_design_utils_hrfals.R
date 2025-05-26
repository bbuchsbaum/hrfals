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
  nb <- nbasis(hrf)
  if (basis_function_index < 1 || basis_function_index > nb) {
    stop("basis_function_index out of range")
  }
  grid <- seq(0, attr(hrf, "span"), by = sampling_frame$TR[1])
  vals <- evaluate(hrf, grid)
  if (is.vector(vals)) {
    vals <- matrix(vals, ncol = 1L)
  }
  kern <- vals[, basis_function_index]
  conv_full <- convolve(raw_timeseries, rev(kern), type = "open")
  conv_full[seq_along(raw_timeseries)]
}

#' Prepare CFALS inputs from a single event term
#'
#' Constructs projected design matrices and related metadata for a
#' specific event term within an `event_model`.
#'
#' @param fmri_data_obj `fmrireg::fmri_dataset` or numeric matrix of BOLD data.
#' @param fmrireg_event_model An `event_model` object.
#' @param target_event_term_name Character name of the event_term to use.
#' @param fmrireg_hrf_basis HRF basis object with `nbasis > 1`.
#' @param confound_obj Optional confound matrix.
#' @param hrf_shape_duration_sec Duration for the HRF reconstruction grid.
#' @param hrf_shape_sample_res_sec Sampling resolution for the HRF grid.
#' @return List with projected design matrices, reconstruction info and
#'   metadata for CFALS engines.
#' @export
prepare_cfals_inputs_from_fmrireg_term <- function(fmri_data_obj,
                                                   fmrireg_event_model,
                                                   target_event_term_name,
                                                   fmrireg_hrf_basis,
                                                   confound_obj = NULL,
                                                   hrf_shape_duration_sec = attr(fmrireg_hrf_basis, "span"),
                                                   hrf_shape_sample_res_sec = fmrireg_event_model$sampling_frame$TR[1]) {
  if (inherits(fmri_data_obj, "fmri_dataset")) {
    Y_raw <- get_data_matrix(fmri_data_obj)
  } else if (is.matrix(fmri_data_obj)) {
    Y_raw <- fmri_data_obj
  } else {
    stop("'fmri_data_obj' must be an 'fmri_dataset' or matrix")
  }

  n_timepoints <- nrow(Y_raw)
  v_voxels <- ncol(Y_raw)
  bad_row_idx <- which(apply(Y_raw, 1, function(r) any(is.na(r))))
  if (length(bad_row_idx) > 0) {
    Y_raw[bad_row_idx, ] <- 0
  }

  all_terms <- fmrireg::terms(fmrireg_event_model)
  if (!target_event_term_name %in% names(all_terms)) {
    stop(sprintf("Target event term '%s' not found in event_model.",
                 target_event_term_name))
  }
  target_term_obj <- all_terms[[target_event_term_name]]
  if (!inherits(target_term_obj, "event_term")) {
    stop("Selected term is not an 'event_term', CF-ALS requires an event_term.")
  }

  X_unconvolved_term <- fmrireg::design_matrix(target_term_obj, drop.empty = TRUE)
  k_conditions <- ncol(X_unconvolved_term)
  if (k_conditions == 0) {
    stop("No estimable conditions in target event term")
  }
  cond_names <- colnames(X_unconvolved_term)

  d_basis_dim <- fmrireg::nbasis(fmrireg_hrf_basis)
  if (d_basis_dim <= 1) {
    stop("CF-ALS requires an hrf_basis with nbasis > 1")
  }

  X_list_raw <- vector("list", k_conditions)
  names(X_list_raw) <- cond_names
  for (c in seq_len(k_conditions)) {
    raw_c <- X_unconvolved_term[, c]
    mat_c <- vapply(seq_len(d_basis_dim), function(j) {
      convolve_timeseries_with_single_basis(raw_c,
                                            fmrireg_hrf_basis,
                                            j,
                                            fmrireg_event_model$sampling_frame)
    }, numeric(n_timepoints))
    mat_c <- matrix(mat_c, nrow = n_timepoints, ncol = d_basis_dim)
    if (length(bad_row_idx) > 0) {
      mat_c[bad_row_idx, ] <- 0
    }
    X_list_raw[[c]] <- mat_c
  }

  time_points_for_shape <- seq(0, hrf_shape_duration_sec, by = hrf_shape_sample_res_sec)
  Phi_recon_matrix <- fmrireg::evaluate(fmrireg_hrf_basis, time_points_for_shape)
  if (is.vector(Phi_recon_matrix) && d_basis_dim == 1) {
    Phi_recon_matrix <- matrix(Phi_recon_matrix, ncol = 1)
  }
  if (ncol(Phi_recon_matrix) != d_basis_dim) {
    stop("Mismatch in Phi_recon_matrix columns and nbasis")
  }
  h_ref_shape_canonical_p_dim <- fmrireg::evaluate(fmrireg::HRF_GLOVER(),
                                                   time_points_for_shape)
  h_ref_shape_canonical_p_dim <- drop(h_ref_shape_canonical_p_dim)
  h_ref_shape_canonical_p_dim <- h_ref_shape_canonical_p_dim / max(abs(h_ref_shape_canonical_p_dim))

  proj <- project_confounds(Y_raw, X_list_raw, confound_obj)
  Y_proj <- proj$Y
  X_list_proj <- proj$X_list

  list(Y_proj = Y_proj,
       X_list_proj = X_list_proj,
       d_basis_dim = d_basis_dim,
       k_conditions = k_conditions,
       Phi_recon_matrix = Phi_recon_matrix,
       h_ref_shape_canonical_p_dim = h_ref_shape_canonical_p_dim,
       n_timepoints = n_timepoints,
       v_voxels = v_voxels,
       bad_row_idx = bad_row_idx,
       target_term_condition_names = cond_names)
}

