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
#' Internal helper that dispatches to [estimate_hrf_cfals()] when a pre-built
#' design list is available. Control parameters are merged with
#' [hrfals_control_defaults()] via [modifyList()].
#'
#' @param y Numeric matrix of BOLD data (time points \eqn{\times} voxels).
#' @param design List produced by [create_cfals_design()] containing at least
#'   `event_model` and `hrf_basis`.
#' @param method Estimation engine to use. Passed to [estimate_hrf_cfals()].
#' @param control List of control parameters overriding defaults.
#' @param ... Additional arguments passed to [estimate_hrf_cfals()].

#' @return An object of class `hrfals_fit`.
#' @export
hrfals_from_design <- function(y, design, method = "cf_als", control = list(), ...) {

  ctrl <- modifyList(hrfals_control_defaults(), control)

  if (is.null(design$event_model) || is.null(design$hrf_basis)) {
    stop("'design' must contain 'event_model' and 'hrf_basis' components")
  }

  target_term <- names(design$event_model$terms)[1]
  penalty_type <- if (is.null(ctrl$R_mat)) "identity" else "custom"
  
  estimate_hrf_cfals(
    fmri_data_obj = y,
    fmrireg_event_model = design$event_model,
    target_event_term_name = target_term,
    hrf_basis_for_cfals = design$hrf_basis,
    method = method,
    lambda_init = ctrl$lambda_init,
    lambda_b = ctrl$lambda_b,
    lambda_h = ctrl$lambda_h,
    lambda_joint = ctrl$lambda_joint,
    R_mat = penalty_type,
    fullXtX = ctrl$fullXtX,
    precompute_xty_flag = ctrl$precompute_xty_flag,
    max_alt = ctrl$max_alt,
    ...
  )
}
