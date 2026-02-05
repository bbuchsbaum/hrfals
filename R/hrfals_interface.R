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
    lambda_s = 0,
    laplacian_obj = NULL,
    h_solver = "auto",
    cg_max_iter = 100,
    cg_tol = 1e-4,
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

  ctrl <- utils::modifyList(hrfals_control_defaults(), control)

  # Check required components
  required_fields <- c("X_list_proj", "Y_proj", "Phi_recon_matrix", 
                      "h_ref_shape_canonical", "d_basis_dim", "k_conditions")
  missing_fields <- setdiff(required_fields, names(design))
  if (length(missing_fields) > 0) {
    stop("'design' must contain the following components: ", 
         paste(missing_fields, collapse = ", "))
  }

  # Extract design components
  Xp <- design$X_list_proj
  Yp <- design$Y_proj  # Use pre-projected data from design
  d <- design$d_basis_dim
  k <- design$k_conditions
  Phi <- design$Phi_recon_matrix
  h_ref_shape_canonical <- design$h_ref_shape_canonical
  n <- design$n_timepoints
  v <- design$v_voxels
  
  # Handle R_mat
  if (is.character(ctrl$R_mat)) {
    R_choice <- match.arg(ctrl$R_mat, c("identity", "basis_default"))
    R_eff <- switch(R_choice,
                    identity = diag(d),
                    basis_default = fmrihrf::penalty_matrix(design$hrf_basis))
  } else if (is.matrix(ctrl$R_mat)) {
    if (nrow(ctrl$R_mat) != d || ncol(ctrl$R_mat) != d) {
      stop(paste("R_mat must be a", d, "x", d, "matrix"))
    }
    R_eff <- ctrl$R_mat
  } else {
    R_eff <- NULL  # Default to NULL to match engine defaults
  }

  # Call the appropriate engine directly
  fit <- switch(method,
    ls_svd_only = ls_svd_engine(Xp, Yp,
                                lambda_init = ctrl$lambda_init,
                                Phi_recon_matrix = Phi,
                                h_ref_shape_canonical = h_ref_shape_canonical),
    ls_svd_1als = ls_svd_1als_engine(Xp, Yp,
                                     lambda_init = ctrl$lambda_init,
                                     lambda_b = ctrl$lambda_b,
                                     lambda_h = ctrl$lambda_h,
                                     lambda_joint = ctrl$lambda_joint,
                                     fullXtX_flag = ctrl$fullXtX,
                                     Phi_recon_matrix = Phi,
                                     h_ref_shape_canonical = h_ref_shape_canonical,
                                     R_mat = R_eff),
	   cf_als = cf_als_engine(Xp, Yp,
	                           lambda_init = ctrl$lambda_init,
	                           lambda_b = ctrl$lambda_b,
	                           lambda_h = ctrl$lambda_h,
	                           lambda_joint = ctrl$lambda_joint,
	                           lambda_s = ctrl$lambda_s,
	                           laplacian_obj = ctrl$laplacian_obj,
                           h_solver = ctrl$h_solver,
                           cg_max_iter = ctrl$cg_max_iter,
                           cg_tol = ctrl$cg_tol,
                           R_mat_eff = R_eff,
                           fullXtX_flag = ctrl$fullXtX,
                           precompute_xty_flag = ctrl$precompute_xty_flag,
                           Phi_recon_matrix = Phi,
                           h_ref_shape_canonical = h_ref_shape_canonical,
                           max_alt = ctrl$max_alt,
                           beta_penalty = list(l1 = 0, alpha = 1, warm_start = TRUE),
                           design_control = list(standardize_predictors = FALSE,
                                                 cache_design_blocks = TRUE))
  )

  # Set condition names if available
	  if (!is.null(design$condition_names)) {
	    rownames(fit$beta) <- design$condition_names
	  }

	  predictor_sds <- design$predictor_sds
	  if (!is.null(predictor_sds) && length(predictor_sds) == nrow(fit$beta)) {
	    if (!is.null(names(predictor_sds)) && !is.null(rownames(fit$beta))) {
	      predictor_sds <- predictor_sds[rownames(fit$beta)]
	    }
	    if (any(!is.finite(predictor_sds)) || any(predictor_sds <= 0)) {
	      stop("design$predictor_sds must be finite and positive when provided")
	    }
	    if (any(abs(predictor_sds - 1) > sqrt(.Machine$double.eps))) {
	      fit$beta <- sweep(fit$beta, 1, predictor_sds, FUN = "/")
	    }
	  }

  # Reconstruct HRF
  recon_hrf <- Phi %*% fit$h

  # Calculate predictions and goodness of fit
  pred_p <- matrix(0, n, v)
  for (c in seq_len(k)) {
    pred_p <- pred_p + (Xp[[c]] %*% fit$h) *
      matrix(rep(fit$beta[c, ], each = n), n, v)
  }
  resids <- Yp - pred_p

	  SST <- colSums((Yp - matrix(colMeans(Yp), n, v, TRUE))^2)
	  SSE <- colSums(resids^2)
	  r2 <- rep(NA_real_, v)
	  ok <- is.finite(SST) & SST > 0
	  r2[ok] <- 1 - SSE[ok] / SST[ok]

  # Create hrfals_fit object
  target_term <- if (!is.null(design$event_model)) {
    names(design$event_model$terms)[1]
  } else {
    "unknown"
  }
  
  hrfals_fit(h_coeffs = fit$h,
             beta_amps = fit$beta,
             method = method,
             lambdas = c(init = ctrl$lambda_init,
                         beta = ctrl$lambda_b,
                         h = ctrl$lambda_h,
                         joint = ctrl$lambda_joint,
                         spatial = ctrl$lambda_s),
             call = match.call(),
             fmrireg_hrf_basis_used = design$hrf_basis,
             target_event_term_name = target_term,
             phi_recon_matrix = Phi,
             design_info = list(d = d, k = k, n = n, v = v, fullXtX = ctrl$fullXtX,
                               predictor_means = if (is.null(design$predictor_means)) rep(0, k) else design$predictor_means,
                               predictor_sds = if (is.null(design$predictor_sds)) rep(1, k) else design$predictor_sds),
             residuals = resids,
             bad_row_idx = if (is.null(design$bad_row_idx)) integer(0) else design$bad_row_idx,
             recon_hrf = recon_hrf,
            gof = r2,
            beta_penalty = list(l1 = 0, alpha = 1, warm_start = TRUE))
}
