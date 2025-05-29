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
#' @param baseline_model Optional baseline model whose design matrix is
#'   projected along with `confound_obj`.
#' @param method Estimation engine to use ("ls_svd_only", "ls_svd_1als", "cf_als").
#' @param lambda_init Ridge penalty for initial LS solve.
#' @param lambda_b Ridge penalty for the beta update. Ignored when
#'   `beta_penalty$l1 > 0` as the Elastic Net solver controls all
#'   regularisation.
#' @param lambda_h Ridge penalty for the h update.
#' @param lambda_joint Joint penalty for the h update.
#' @param lambda_s Spatial regularization strength controlling the amount of
#'   voxel-wise smoothing.  The effective physical smoothness depends on the
#'   value of `lambda_s` relative to the voxel size.  For similar behaviour
#'   across datasets with different resolutions consider scaling
#'   `lambda_s` by `1 / mean(voxel_size)^2`.  Set to zero to disable spatial
#'   smoothing.
#' @param laplacian_obj Optional list containing the Laplacian matrix `L` and
#'   degree vector `degree` as returned by [build_voxel_laplacian()].
#'   Required when `lambda_s > 0`.
#' @param h_solver Solver to use for the spatial h-update. One of
#'   "direct", "cg" or "auto".
#' @param cg_max_iter Maximum iterations for the conjugate gradient solver when
#'   `h_solver = "cg"`.
#' @param cg_tol Convergence tolerance for the conjugate gradient solver.
#' @param R_mat Either the character string "identity" (default) or
#'   "basis_default" to indicate how the penalty matrix should be
#'   constructed, or a numeric matrix providing a custom penalty for the
#'   h update.
#' @param fullXtX Logical. If `TRUE`, the h-update step uses the full
#'   Gramian matrix \eqn{(\sum_l \beta_l X_l)^\top (\sum_m \beta_m X_m)}
#'   with cross-condition terms. If `FALSE` (default), the Gramian is
#'   approximated by omitting cross-condition terms,
#'   \eqn{\sum_l \beta_l^2 X_l^\top X_l}. A single shared HRF
#'   coefficient vector is still estimated per voxel.
#' @param precompute_xty_flag Logical; passed to `cf_als_engine`.
#' @param max_alt Number of alternating updates for `cf_als`.
#' @param beta_penalty List with elements `l1`, `alpha`, and `warm_start`
#'   passed to `cf_als_engine` for sparse beta estimation. Set `l1 > 0`
#'   to enable an Elastic Net update.
#' @param design_control List of design matrix processing options. Set
#'   `standardize_predictors = TRUE` to z-score continuous predictors
#'   before estimation.
#' @param hrf_shape_duration Duration in seconds for reconstructed HRF grid.
#' @param hrf_shape_resolution Sampling resolution of the HRF grid.
#' @return An `hrfals_fit` object.
#' @export
estimate_hrf_cfals <- function(fmri_data_obj,
                               fmrireg_event_model,
                               target_event_term_name,
                               hrf_basis_for_cfals,
                               confound_obj = NULL,
                               baseline_model = NULL,
                               method = c("ls_svd_1als", "ls_svd_only", "cf_als"),
                               lambda_init = 1,
                               lambda_b = 10,
                               lambda_h = 1,
                               lambda_joint = 0,
                               lambda_s = 0,
                               laplacian_obj = NULL,
                               h_solver = c("direct", "cg", "auto"),
                               cg_max_iter = 100,
                               cg_tol = 1e-4,
                               R_mat = c("identity", "basis_default"),
                               fullXtX = FALSE,
                               precompute_xty_flag = TRUE,
                               max_alt = 10,
                               beta_penalty = list(l1 = 0, alpha = 1, warm_start = TRUE),
                               design_control = list(standardize_predictors = TRUE,
                                                     cache_design_blocks = TRUE),
                               hrf_shape_duration = attr(hrf_basis_for_cfals, "span"),
                               hrf_shape_resolution = fmrireg_event_model$sampling_frame$TR[1],
                               ...) {
  method <- match.arg(method)

  prep <- create_cfals_design(fmri_data_obj,
                             fmrireg_event_model,
                             hrf_basis_for_cfals,
                             confound_obj = confound_obj,
                             baseline_model = baseline_model,
                             hrf_shape_duration_sec = hrf_shape_duration,
                             hrf_shape_sample_res_sec = hrf_shape_resolution,
                             design_control = design_control)

  Xp <- prep$X_list_proj
  Yp <- prep$Y_proj
  d <- prep$d_basis_dim
  k <- prep$k_conditions
  Phi <- prep$Phi_recon_matrix
  h_ref_shape_canonical <- prep$h_ref_shape_canonical
  n <- prep$n_timepoints
  v <- prep$v_voxels

  if (is.character(R_mat)) {
    R_choice <- match.arg(R_mat, c("identity", "basis_default"))
    R_eff <- switch(R_choice,
                    identity = diag(d),
                    basis_default = fmrireg::penalty_matrix(hrf_basis_for_cfals))
  } else if (is.matrix(R_mat)) {
    if (nrow(R_mat) != d || ncol(R_mat) != d) {
      stop(paste("R_mat must be a", d, "x", d, "matrix"))
    }
    R_eff <- R_mat
  } else {
    stop("R_mat must be 'identity', 'basis_default', or a numeric matrix")
  }

  fit <- switch(method,
    ls_svd_only = ls_svd_engine(Xp, Yp,
                                lambda_init = lambda_init,
                                Phi_recon_matrix = Phi,
                                h_ref_shape_canonical = h_ref_shape_canonical),
    ls_svd_1als = ls_svd_1als_engine(Xp, Yp,
                                     lambda_init = lambda_init,
                                     lambda_b = lambda_b,
                                     lambda_h = lambda_h,
                                     lambda_joint = lambda_joint,
                                     fullXtX_flag = fullXtX,
                                     Phi_recon_matrix = Phi,
                                     h_ref_shape_canonical = h_ref_shape_canonical,
                                     R_mat = R_eff),
   cf_als = cf_als_engine(Xp, Yp,
                           lambda_b = lambda_b,
                           lambda_h = lambda_h,
                           lambda_joint = lambda_joint,
                           lambda_s = lambda_s,
                           laplacian_obj = laplacian_obj,
                           h_solver = h_solver,
                           cg_max_iter = cg_max_iter,
                           cg_tol = cg_tol,
                           R_mat_eff = R_eff,
                           fullXtX_flag = fullXtX,
                           precompute_xty_flag = precompute_xty_flag,
                           Phi_recon_matrix = Phi,
                           h_ref_shape_canonical = h_ref_shape_canonical,
                           max_alt = max_alt,
                           beta_penalty = beta_penalty,
                           design_control = design_control)
  )

  rownames(fit$beta) <- prep$condition_names

  if (isTRUE(design_control$standardize_predictors)) {
    fit$beta <- sweep(fit$beta, 1, prep$predictor_sds, FUN = "/")
  }
  recon_hrf <- Phi %*% fit$h

  pred_p <- matrix(0, n, v)
  for (c in seq_len(k)) {
    pred_p <- pred_p + (Xp[[c]] %*% fit$h) *
      matrix(rep(fit$beta[c, ], each = n), n, v)
  }
  resids <- Yp - pred_p

  SST <- colSums((Yp - matrix(colMeans(Yp), n, v, TRUE))^2)
  SSE <- colSums(resids^2)
  r2 <- 1 - SSE / SST

  hrfals_fit(h_coeffs = fit$h,
             beta_amps = fit$beta,
             method = method,
             lambdas = c(init = lambda_init,
                         beta = lambda_b,
                         h = lambda_h,
                         joint = lambda_joint,
                         spatial = lambda_s),
             call = match.call(),
             fmrireg_hrf_basis_used = hrf_basis_for_cfals,
             target_event_term_name = target_event_term_name,
             phi_recon_matrix = Phi,
             design_info = list(d = d, k = k, n = n, v = v, fullXtX = fullXtX,
                               predictor_means = prep$predictor_means,
                               predictor_sds = prep$predictor_sds),
             residuals = resids,
             bad_row_idx = prep$bad_row_idx,
             recon_hrf = recon_hrf,
            gof = r2)
}

#' Spatially-Regularised CF-ALS Convenience Wrapper
#'
#' Calls [estimate_hrf_cfals()] with a predefined non-zero `lambda_s`.
#'
#' @inheritParams estimate_hrf_cfals
#' @param lambda_s_default Default spatial regularisation strength.
#' @export
estimate_hrf_spatial_cfals <- function(fmri_data_obj,
                                       fmrireg_event_model,
                                       target_event_term_name,
                                       hrf_basis_for_cfals,
                                       laplacian_obj,
                                       lambda_s_default = 0.1,
                                       ...) {
  if (missing(laplacian_obj)) {
    stop("laplacian_obj with elements L and degree must be provided")
  }
  
  estimate_hrf_cfals(fmri_data_obj,
                     fmrireg_event_model,
                     target_event_term_name,
                     hrf_basis_for_cfals,
                     laplacian_obj = laplacian_obj,
                     lambda_s = lambda_s_default,
                     ...)
}

