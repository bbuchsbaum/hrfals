#' Beta-series estimation using fast LSS
#'
#' Builds trial regressors from a CF-ALS fit and event information and
#' calls \code{fastlss_shared()} or \code{fastlss_voxel()} as appropriate.
#'
#' @param cf_fit An object of class `hrfals_fit` as returned by
#'   [hrfals()] or [fmrireg_cfals()].
#' @param events A `data.frame` with onset information or an
#'   `fmrireg::event_model`.
#' @param fmri_data_obj BOLD data matrix or `fmrireg::fmri_dataset`
#'   matching the data used for CF-ALS.
#' @param confound_obj Optional confound matrix with the same number of
#'   rows as `fmri_data_obj`.
#' @param baseline_model Optional baseline model for confound projection.
#' @param mode Character string specifying the LSS kernel variant.
#'   If "auto" (default) the function selects "shared" when
#'   a single HRF is present in `cf_fit$h_coeffs` and "voxel" otherwise.
#' @param whitening_matrix Optional whitening matrix applied to all
#'   design matrices before running the kernel.
#' @param ... Additional arguments passed to [fastlss_shared()] or
#'   [fastlss_voxel()].
#' @return An object of class `fastlss_fit` containing the trial
#'   coefficient matrix and metadata.
#' @export
hrfals_lss <- function(cf_fit, events, fmri_data_obj,
                       confound_obj = NULL,
                       baseline_model = NULL,
                       mode = c("auto", "shared", "voxel"),
                       whitening_matrix = NULL,
                       ...) {
  if (!inherits(cf_fit, "hrfals_fit"))
    stop("'cf_fit' must be an object of class 'hrfals_fit'")

  mode <- match.arg(mode)

  if (inherits(events, "event_model")) {
    event_model <- events
  } else if (is.data.frame(events)) {
    sframe <- if (inherits(fmri_data_obj, "fmri_dataset"))
      fmri_data_obj$sampling_frame
    else
      fmrireg::sampling_frame(nrow(fmri_data_obj),
                              TR = attr(fmri_data_obj, "TR")[1])
    block_formula <- if ("block" %in% names(events)) ~ block else NULL
    event_model <- fmrireg::event_model(onset ~ hrf(condition),
                                        data = events,
                                        block = block_formula,
                                        sampling_frame = sframe)
  } else {
    stop("'events' must be a data.frame or an 'event_model'")
  }

  basis <- cf_fit$fmrireg_hrf_basis_used
  if (is.null(basis))
    stop("cf_fit must contain 'fmrireg_hrf_basis_used'")

  design <- create_cfals_design(fmri_data_obj, event_model, basis,
                                confound_obj = confound_obj,
                                baseline_model = baseline_model)
  Y <- design$Y_proj
  X_list <- design$X_list_proj

  if (mode == "auto") {
    mode <- if (ncol(cf_fit$h_coeffs) > 1) "voxel" else "shared"
  }

  if (mode == "shared") {
    h_shared <- if (ncol(cf_fit$h_coeffs) > 1)
      rowMeans(cf_fit$h_coeffs) else drop(cf_fit$h_coeffs[, 1])
    C <- matrix(0, nrow(Y), length(X_list))
    for (t in seq_along(X_list)) {
      C[, t] <- X_list[[t]] %*% h_shared
    }
    p_vec <- rep(0, nrow(Y))
    B <- fastlss_shared(Y, matrix(0, nrow(Y), 0), C, p_vec,
                    W = whitening_matrix, ...)
  } else {
    p_vec <- rep(0, nrow(Y))
    B <- fastlss_voxel(Y, matrix(0, nrow(Y), 0), X_list,
                    cf_fit$h_coeffs, p_vec,
                    W = whitening_matrix, ...)
  }

  dimnames(B) <- list(names(X_list), colnames(Y))
  fastlss_fit(B, mode = mode, cfals_fit = cf_fit,
              events = events, hrf_basis = basis,
              call = match.call(),
              whitening_matrix = whitening_matrix)
}
