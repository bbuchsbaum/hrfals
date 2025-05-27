#' hrfals control defaults
#'
#' Returns a list of default tuning parameters for `hrfals()`.
#'
#' @return Named list of defaults used by [hrfals()].
#' @export
hrfals_control_defaults <- function() {
  list(
    lambda_init = 1,
    lambda_b = 10,
    lambda_h = 1,
    lambda_joint = 0,
    R_mat = NULL,
    fullXtX = FALSE,
    precompute_xty_flag = TRUE,
    max_alt = 1
  )
}

#' Fit HRFs using CF-ALS (Design Interface)
#'
#' Internal helper that dispatches to [fmrireg_cfals()] when a pre-built
#' design list is available. Control parameters are merged with
#' [hrfals_control_defaults()] via [modifyList()].
#'
#' @param y Numeric matrix of BOLD data (time points \eqn{\times} voxels).
#' @param design List produced by [create_cfals_design()] containing at least
#'   `event_model` and `hrf_basis`.
#' @param method Estimation engine to use. Passed to [fmrireg_cfals()].
#' @param control List of control parameters overriding defaults.
#' @param ... Additional arguments passed to [fmrireg_cfals()].
#' @return An object of class `fmrireg_cfals_fit`.
#' @noRd
hrfals_from_design <- function(y, design, method = "cf_als", control = list(), ...) {
  ctrl <- modifyList(hrfals_control_defaults(), control)

  if (is.null(design$event_model) || is.null(design$hrf_basis)) {
    stop("'design' must contain 'event_model' and 'hrf_basis' components")
  }

  fmrireg_cfals(
    fmri_data_obj = y,
    event_model = design$event_model,
    hrf_basis = design$hrf_basis,
    method = method,
    lambda_init = ctrl$lambda_init,
    lambda_b = ctrl$lambda_b,
    lambda_h = ctrl$lambda_h,
    lambda_joint = ctrl$lambda_joint,
    R_mat = ctrl$R_mat,
    fullXtX = ctrl$fullXtX,
    precompute_xty_flag = ctrl$precompute_xty_flag,
    max_alt = ctrl$max_alt,
    ...
  )
}
