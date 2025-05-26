#' Estimate HRF for a target event term using CF-ALS
#'
#' High level wrapper around the CFALS engines operating on a single
#' `event_term` within an `fmrireg::event_model`. Design matrices and
#' projection are handled by `create_cfals_design`.
#'
#' @param fmri_data_obj `fmrireg::fmri_dataset` or numeric BOLD matrix.
#' @param fmrireg_event_model An `event_model` describing the full design.
#' @param target_event_term_name Name of the event_term to estimate.
#' @param hrf_basis_for_cfals An `HRF` object with `nbasis > 1`.
#' @param confound_obj Optional confound matrix.
#' @param method Estimation engine to use ("ls_svd_only", "ls_svd_1als", "cf_als").
#' @param lambda_init Ridge penalty for initial LS solve.
#' @param lambda_b Ridge penalty for the beta update.
#' @param lambda_h Ridge penalty for the h update.
#' @param penalty_R_mat_type How to construct the penalty matrix. One of
#'   "identity", "basis", or "custom". If "custom", supply `R_mat`.
#' @param R_mat Optional custom penalty matrix for the h update.
#' @param fullXtX Logical; passed to the estimation engine.
#' @param precompute_xty_flag Logical; passed to `cf_als_engine`.
#' @param max_alt Number of alternating updates for `cf_als`.
#' @param hrf_shape_duration Duration in seconds for reconstructed HRF grid.
#' @param hrf_shape_resolution Sampling resolution of the HRF grid.
#' @return An `hrfals_fit` object.
#' @export
estimate_hrf_cfals <- function(fmri_data_obj,
                               fmrireg_event_model,
                               target_event_term_name,
                               hrf_basis_for_cfals,
                               confound_obj = NULL,
                               method = c("ls_svd_1als", "ls_svd_only", "cf_als"),
                               lambda_init = 1,
                               lambda_b = 10,
                               lambda_h = 1,
                               penalty_R_mat_type = c("identity", "basis", "custom"),
                               R_mat = NULL,
                               fullXtX = FALSE,
                               precompute_xty_flag = TRUE,
                               max_alt = 1,
                               hrf_shape_duration = attr(hrf_basis_for_cfals, "span"),
                               hrf_shape_resolution = fmrireg_event_model$sampling_frame$TR[1],
                               ...) {
  method <- match.arg(method)
  penalty_R_mat_type <- match.arg(penalty_R_mat_type)

  prep <- create_cfals_design(fmri_data_obj,
                             fmrireg_event_model,
                             hrf_basis_for_cfals,
                             confound_obj = confound_obj,
                             hrf_shape_duration_sec = hrf_shape_duration,
                             hrf_shape_sample_res_sec = hrf_shape_resolution)

  Xp <- prep$X_list_proj
  Yp <- prep$Y_proj
  d <- prep$d_basis_dim
  k <- prep$k_conditions
  Phi <- prep$Phi_recon_matrix
  n <- prep$n_timepoints
  v <- prep$v_voxels

  R_eff <- switch(penalty_R_mat_type,
                  identity = diag(d),
                  basis = penalty_matrix(hrf_basis_for_cfals),
                  custom = {
                    if (is.null(R_mat)) stop("R_mat must be supplied for custom penalty")
                    R_mat
                  })

  fit <- switch(method,
    ls_svd_only = ls_svd_engine(Xp, Yp,
                                lambda_init = lambda_init,
                                h_ref_shape_norm = NULL,
                                R_mat = R_eff),
    ls_svd_1als = ls_svd_1als_engine(Xp, Yp,
                                     lambda_init = lambda_init,
                                     lambda_b = lambda_b,
                                     lambda_h = lambda_h,
                                     fullXtX_flag = fullXtX,
                                     h_ref_shape_norm = NULL,
                                     R_mat = R_eff),
    cf_als = cf_als_engine(Xp, Yp,
                           lambda_b = lambda_b,
                           lambda_h = lambda_h,
                           R_mat_eff = R_eff,
                           fullXtX_flag = fullXtX,
                           precompute_xty_flag = precompute_xty_flag,
                           h_ref_shape_norm = NULL,
                           max_alt = max_alt)
  )

  rownames(fit$beta) <- prep$condition_names
  recon_hrf <- Phi %*% fit$h

  pred_p <- Reduce(`+`, Map(function(Xc, bc) {
    Xc %*% (fit$h * matrix(bc, nrow = d, ncol = v, byrow = TRUE))
  }, Xp, asplit(fit$beta, 1)))
  resids <- Yp - pred_p

  SST <- colSums((Yp - matrix(colMeans(Yp), n, v, TRUE))^2)
  SSE <- colSums(resids^2)
  r2 <- 1 - SSE / SST

  hrfals_fit(h_coeffs = fit$h,
             beta_amps = fit$beta,
             method = method,
             lambdas = c(init = lambda_init,
                         beta = lambda_b,
                         h = lambda_h),
             call = match.call(),
             fmrireg_hrf_basis_used = hrf_basis_for_cfals,
             target_event_term_name = target_event_term_name,
             phi_recon_matrix = Phi,
             design_info = list(d = d, k = k, n = n, v = v, fullXtX = fullXtX),
             residuals = resids,
             bad_row_idx = prep$bad_row_idx,
             recon_hrf = recon_hrf,
             gof = r2)
}

