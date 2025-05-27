#' Create CFALS design from events and data
#'
#' Convenience wrapper that first constructs an `fmrireg::event_model` from a
#' data.frame of events and then calls [create_cfals_design()] to generate
#' the matrices required for CFALS estimation.
#'
#' @param events Data frame with at least an `onset` column and any factors used
#'   in the event model. If a `block` column is not present a single block is
#'   assumed.
#' @param TR Numeric. Repetition time in seconds.
#' @param basis An `HRF` basis object.
#' @param fmri_data_obj BOLD data matrix (time points \eqn{\times} voxels) or
#'   `fmrireg::fmri_dataset`.
#' @param confound_obj Optional confound matrix with the same number of rows as
#'   `fmri_data_obj`.
#' @param baseline_model Optional baseline model whose design matrix is
#'   projected alongside `confound_obj`.
#' @param formula Model formula passed to [fmrireg::event_model()]. The default
#'   expects a column named `condition` in `events`.
#' @param block Formula specifying the block column. Defaults to `~ block` if a
#'   `block` column is present.
#' @param ... Additional arguments passed to [fmrireg::event_model()].
#'
#' @return A list as returned by [create_cfals_design()].
#' @export
hrfals_design <- function(events, TR, basis, fmri_data_obj,
                          confound_obj = NULL,
                          baseline_model = NULL,
                          formula = onset ~ hrf(condition),
                          block = if ("block" %in% names(events)) ~ block else NULL,
                          ...) {
  if (inherits(fmri_data_obj, "fmri_dataset")) {
    sframe <- fmri_data_obj$sampling_frame
  } else if (is.matrix(fmri_data_obj)) {
    sframe <- fmrireg::sampling_frame(nrow(fmri_data_obj), TR = TR)
  } else {
    stop("'fmri_data_obj' must be an 'fmri_dataset' or matrix")
  }

  emod <- fmrireg::event_model(formula = formula,
                               data = events,
                               block = block,
                               sampling_frame = sframe,
                               ...)

  create_cfals_design(fmri_data_obj,
                      emod,
                      basis,
                      confound_obj = confound_obj,
                      baseline_model = baseline_model)
}
