#' Estimate Rank-1 HRF Using LS+SVD/CF-ALS Methods
#'
#' High level wrapper that prepares design matrices and dispatches to
#' the desired estimation engine.  This function supports the fast
#' \emph{LS+SVD} initialisation, the one-step refinement
#' \emph{LS+SVD+1ALS}, or the full alternating scheme implemented in
#' `cf_als_engine`.
#'
#' @param fmri_data_obj An `fmri_dataset` or numeric matrix of BOLD data
#'   (time points x voxels). If a dataset, sampling information is
#'   taken from the object.
#' @param event_model An `event_model` describing the stimuli to use
#'   for HRF estimation.
#' @param hrf_basis An `HRF` basis object used for the convolution
#'   design matrices.
#' @param confound_obj Optional matrix of confound regressors with the
#'   same number of rows as the data matrix.
#' @param method Estimation method. One of "ls_svd_only",
#'   "ls_svd_1als" (default) or "cf_als".
#' @param lambda_init Ridge penalty for the initial GLM solve used by
#'   `ls_svd` based methods.
#' @param lambda_b Ridge penalty for the beta update step. This value
#'   is ignored when `beta_penalty$l1 > 0` because the Elastic Net
#'   solver handles all regularisation.
#' @param lambda_h Ridge penalty for the h update step.
#' @param lambda_joint Joint ridge penalty applied to both beta and h updates.
#'   This helps prevent see-saw effects between the two parameter blocks.
#'   Recommended range: 0.2-2. Default is 0 (disabled).
#' @param R_mat Optional penalty matrix for the h (HRF coefficient) update.
#'   If `NULL` (default), an identity matrix is used, corresponding to a simple
#'   ridge penalty. If a basis-specific penalty (e.g., for smoothing) is
#'   available from `hrf_basis`, it can be passed here.
#' @param fullXtX Logical. If `TRUE`, the h-update step uses the full
#'   Gramian matrix \eqn{(\sum_l \beta_l X_l)^\top (\sum_m \beta_m X_m)},
#'   including cross-condition terms \eqn{\beta_l \beta_m X_l^\top X_m}.
#'   If `FALSE` (default), cross-condition terms are omitted and the
#'   Gramian is approximated by \eqn{\sum_l \beta_l^2 X_l^\top X_l}.
#'   In both cases a single shared HRF coefficient vector is estimated
#'   per voxel.
#' @param precompute_xty_flag Logical; passed to `cf_als_engine` to control
#'   precomputation of `XtY` matrices.
#' @param max_alt Number of alternating updates after initialisation
#'   when `method = "cf_als"`.
#' @param ... Additional arguments passed to the underlying estimation engine.
#' @return An object of class `hrfals_fit` containing the
#'   estimated HRF coefficients and amplitudes.
#' @details
#' The `method` argument selects between the closed-form
#' \code{"ls_svd_only"}, the default \code{"ls_svd_1als"} which adds one
#' ALS refinement step, or the iterative \code{"cf_als"} engine.  The
#' ridge penalties \code{lambda_init}, \code{lambda_b} and
#' \code{lambda_h} control regularisation of the initial solve, the
#' beta-update and the h-update respectively.  When
#' \code{fullXtX = TRUE} the h-update uses the full cross-condition
#' Gramian; otherwise an approximation ignoring off-diagonal terms is
#' used. R\eqn{^2} is computed on the
#' data after confound projection.
#'
#' @examples
#' \dontrun{
#' library(fmrireg)
#' 
#' # Create sampling frame and event model
#' sframe <- fmrireg::sampling_frame(blocklens = 40, TR = 1)
#' ev_df <- data.frame(onset = c(5, 15, 25), block = 1, cond = "A")
#' emod <- fmrireg::event_model(onset ~ hrf(cond), data = ev_df, 
#'                              block = ~ block, sampling_frame = sframe)
#' 
#' # Simulate some BOLD data
#' Y_matrix <- matrix(rnorm(40 * 5), 40, 5) # 40 timepoints, 5 voxels
#' 
#' # Fit using CF-ALS with SPMG3 basis (3 basis functions)
#' fit <- fmrireg_cfals(Y, emod, HRF_SPMG3)
#' print(fit)
#' }
#' @keywords internal
fmrireg_cfals <- function(fmri_data_obj,
                         event_model,
                         hrf_basis,
                         confound_obj = NULL,
                         method = c("ls_svd_1als", "ls_svd_only", "cf_als"),
                         lambda_init = 1,
                         lambda_b = 10,
                         lambda_h = 1,
                         lambda_joint = 0,
                         R_mat = NULL,
                         fullXtX = FALSE,
                         precompute_xty_flag = TRUE,
                         max_alt = 1,
                         ...) {

  method <- match.arg(method)

  if (inherits(fmri_data_obj, "fmri_dataset")) {
    Y <- get_data_matrix(fmri_data_obj)
  } else if (is.matrix(fmri_data_obj)) {
    Y <- fmri_data_obj
  } else {
    stop("'fmri_data_obj' must be an 'fmri_dataset' or matrix")
  }

  sframe <- if (inherits(fmri_data_obj, "fmri_dataset"))
    fmri_data_obj$sampling_frame else attr(fmri_data_obj, "sampling_frame")

  if (is.null(sframe)) {
    stop("Sampling information could not be determined from input")
  }

  design <- create_cfals_design(fmri_data_obj,
                                event_model,
                                hrf_basis,
                                confound_obj = confound_obj)

  Xp <- design$X_list_proj
  Yp <- design$Y_proj
  cond_names <- design$condition_names

  fit <- switch(method,
    ls_svd_only = ls_svd_engine(Xp, Yp,
                                lambda_init = lambda_init,
                                Phi_recon_matrix = design$Phi_recon_matrix,
                                h_ref_shape_canonical = design$h_ref_shape_canonical,
                                ...),
    ls_svd_1als = ls_svd_1als_engine(Xp, Yp,
                                     lambda_init = lambda_init,
                                     lambda_b = lambda_b,
                                     lambda_h = lambda_h,
                                     lambda_joint = lambda_joint,
                                     fullXtX_flag = fullXtX,
                                     Phi_recon_matrix = design$Phi_recon_matrix,
                                     h_ref_shape_canonical = design$h_ref_shape_canonical,
                                     R_mat = R_mat, ...),
    cf_als = cf_als_engine(Xp, Yp,
                           lambda_b = lambda_b,
                           lambda_h = lambda_h,
                           lambda_joint = lambda_joint,
                           R_mat_eff = R_mat,
                           fullXtX_flag = fullXtX,
                           precompute_xty_flag = precompute_xty_flag,
                           Phi_recon_matrix = design$Phi_recon_matrix,
                           h_ref_shape_canonical = design$h_ref_shape_canonical,
                           max_alt = max_alt, ...)
  )

  rownames(fit$beta) <- cond_names

  Phi <- design$Phi_recon_matrix
  recon_hrf <- Phi %*% fit$h

  n <- nrow(Yp)
  v <- ncol(Yp)
  k <- length(Xp)
  d <- nrow(fit$h)
  pred_p <- matrix(0, n, v)
  for (c in seq_len(k)) {
    pred_p <- pred_p + (Xp[[c]] %*% fit$h) *
      matrix(rep(fit$beta[c, ], each = n), n, v)
  }
  resids <- Yp - pred_p

  SST <- colSums((Yp - matrix(colMeans(Yp), n, v, TRUE))^2)
  SSE <- colSums(resids^2)
  r2 <- 1 - SSE / SST

  out <- hrfals_fit(h_coeffs = fit$h,
                           beta_amps = fit$beta,
                           method = method,
                           lambdas = c(init = lambda_init,
                                       beta = lambda_b,
                                       h = lambda_h,
                                       joint = lambda_joint),
                           call = match.call(),
                           fmrireg_hrf_basis_used = hrf_basis,
                           design_info = list(d = design$d,
                                              k = design$k,
                                              n = n,
                                              v = v,
                                              fullXtX = fullXtX),
                           residuals = resids,
                           recon_hrf = recon_hrf,
                           gof = r2)
  out
}

#' Fit Rank-1 HRF Using CF-ALS
#'
#' Convenience wrapper for the original CF-ALS implementation.  This
#' function simply calls [fmrireg_cfals()] with `method = "cf_als"` and
#' retains the historical argument names.
#'
#' @param fmri_data_obj An `fmri_dataset` or numeric matrix of BOLD data
#'   (time points x voxels). If a dataset, sampling information is
#'   taken from the object.
#' @param event_model An `event_model` describing the stimuli to use
#'   for HRF estimation.
#' @param hrf_basis An `HRF` basis object from the `fmrireg` package
#'   (e.g., `HRF_BSPLINE`, `HRF_TENT`, `HRF_SPMG3`). This determines the
#'   set of basis functions used to model the HRF. The CF-ALS method
#'   can work with any basis where `nbasis > 1`.
#' @param confound_obj Optional numeric matrix of confound regressors (time
#'   points x number of confounds). These nuisance variables are projected
#'   out from the BOLD data and design matrices prior to CF-ALS estimation.
#'   Defaults to `NULL` (no confounds).
#' @param lam_beta Numeric. The regularization parameter for the beta (amplitude)
#'   update step. Controls the L2 penalty on the amplitude coefficients.
#'   Defaults to 10.
#' @param lam_h Numeric. The regularization parameter for the h (HRF coefficient)
#'   update step. Controls the L2 penalty on the HRF basis coefficients.
#'   Defaults to 1.
#' @param R_mat Optional penalty matrix for the h (HRF coefficient) update.
#'   If `NULL` (default), an identity matrix is used, corresponding to a simple
#'   ridge penalty. If a basis-specific penalty (e.g., for smoothing) is
#'   available from `hrf_basis`, it can be passed here.
#' @param fullXtX Logical. If `TRUE`, the h-update step uses the full
#'   Gramian matrix \eqn{(\sum_l \beta_l X_l)^\top (\sum_m \beta_m X_m)},
#'   including cross-condition terms \eqn{\beta_l \beta_m X_l^\top X_m}.
#'   If `FALSE` (default), cross-condition terms are omitted and the
#'   Gramian is approximated by \eqn{\sum_l \beta_l^2 X_l^\top X_l}. In
#'   both cases a single shared HRF coefficient vector is estimated per
#'   voxel.
#' @param max_alt Integer. The maximum number of alternations between the beta
#'   and h updates. The proposal notes that empirically one alternation
#'   (`max_alt = 1`) after SVD initialization is often sufficient.
#'   Defaults to 1.
#' @param ... Additional arguments passed to the underlying estimation engine.
#'
#' @return An object of class `hrfals_fit`. This object contains:
#'   \itemize{
#'     \item `h`: A d x v matrix of HRF basis coefficients (d = number of basis functions, v = number of voxels).
#'     \item `beta`: A k x v matrix of condition amplitudes (k = number of conditions).
#'     \item `reconstructed_hrfs`: A p x v matrix of the actual HRF shapes (p = length of HRF, v = number of voxels), reconstructed using the `hrf_basis` and the estimated `h` coefficients.
#'     \item `residuals`: An n x v matrix of model residuals after fitting (n = number of timepoints).
#'     \item `hrf_basis_used`: The `hrf_basis` object that was supplied to the function.
#'     \item `lambdas_used`: A named list or vector containing the regularization parameters (`lam_beta`, `lam_h`) used in the estimation.
#'     \item `design_info`: A list containing dimensions and flags used during estimation (e.g., d, k, n, v, `fullXtX`).
#'     \item `gof_per_voxel`: Optional. Goodness-of-fit statistics per voxel, such as R-squared.
#'     \item (Other elements as defined by the `hrfals_fit` class structure from the proposal, like `call`).
#'   }
#'
#' @details
#' This function implements the Confound-Free Alternating Least Squares (CF-ALS)
#' algorithm for data-driven estimation of Hemodynamic Response Functions (HRFs)
#' from fMRI data. It estimates HRF coefficients (`h`) and activation
#' amplitudes (`beta`) simultaneously using a rank-1 decomposition model:
#' Y ≈ D(h)β^T, where D(h) is the design matrix formed by convolving
#' stimulus onsets with the HRF (`h`).
#'
#' Key steps of the algorithm include:
#' \enumerate{
#'   \item Confound Projection: Nuisance regressors (if provided via `confound_obj`) are removed from both the BOLD data (Y) and the HRF basis design matrices (X.list) using QR-based orthogonal projection.
#'   \item Precomputation: Quantities like X^T X and X^T Y are precomputed for efficiency.
#'   \item SVD Initialization: Robust starting values for `h` and `beta` are obtained using a regularized least squares solution followed by Singular Value Decomposition (SVD) of the initial coefficient estimates. (This is handled by the main `fmrireg_cfals` function when `method="cf_als"`).
#'   \item CF-ALS Alternation: The algorithm alternates between updating `beta` (amplitudes) holding `h` (HRF coefficients) fixed, and updating `h` holding `beta` fixed. This process is repeated for `max_alt` iterations.
#'     \itemize{
#'       \item β-update: Solves (G_mat + λ_β I)β_v = D(h)^T y_v for each voxel v, where G_mat = D(h)^T D(h).
#'       \item h-update: Solves (LHS + λ_h I)h = RHS, where LHS and RHS depend on the data and current beta estimates. The `fullXtX` parameter influences the construction of LHS.
#'     }
#'   \item Identifiability: The estimated HRF coefficients `h` are typically normalized (e.g., ||Φh||_∞ = 1, where Φ is the basis reconstruction matrix) and their sign aligned with a canonical HRF shape to ensure consistent scaling and polarity across voxels/conditions.
#' }
#'
#' The function is designed to be compatible with any HRF basis defined in the
#' `fmrireg` package (provided `nbasis(hrf_basis) > 1`). The `R_mat`
#' parameter allows for basis-specific penalty matrices for HRF coefficient
#' regularization, although the default is an identity matrix (standard ridge penalty).
#'
#' R\eqn{^2} (coefficient of determination) is computed on the data *after*
#' any confound projection has been applied.
#'
#' @seealso [fmrireg_cfals()] for the more general wrapper allowing different estimation methods.
#' @references (If applicable, add references to papers describing the CF-ALS method or its implementation).
#'
#' @examples
#' # Simulate data
#' sframe <- fmrireg::sampling_frame(blocklens = 40, TR = 1)
#' ev_df <- data.frame(onset = c(5, 15, 25), block = 1, cond = "A")
#' emod <- fmrireg::event_model(onset ~ hrf(cond, basis = fmrireg::HRF_SPMG3),
#'                              data = ev_df, block = ~ block,
#'                              sampling_frame = sframe)
#' Y_matrix <- matrix(rnorm(40 * 5), 40, 5) # 40 timepoints, 5 voxels
#'
#' # Fit using fmrireg_hrf_cfals
#' cfals_fit <- fmrireg_hrf_cfals(
#'   fmri_data_obj = Y_matrix,
#'   event_model = emod,
#'   hrf_basis = fmrireg::HRF_SPMG3, # Using SPMG3 basis (3 basis functions)
#'   lam_beta = 5,
#'   lam_h = 0.5,
#'   max_alt = 1
#' )
#'
#' print(cfals_fit)
#' summary(cfals_fit) # If a summary method is defined
#' # plot(cfals_fit)    # If a plot method is defined
#'
#' # Example with B-spline basis
#' bspline_basis <- fmrireg::HRF_BSPLINE(knots = c(0, 4, 8, 12, 20))
#' emod_bspline <- fmrireg::event_model(onset ~ hrf(cond, basis = bspline_basis),
#'                                      data = ev_df, block = ~ block,
#'                                      sampling_frame = sframe)
#' cfals_fit_bspline <- fmrireg_hrf_cfals(
#'   fmri_data_obj = Y_matrix,
#'   event_model = emod_bspline,
#'   hrf_basis = bspline_basis,
#'   max_alt = 2
#' )
#' print(cfals_fit_bspline)
#'
#' Estimate HRF Using CF-ALS
#'
#' Primary user-facing wrapper for Confound-Free Alternating Least Squares.
#' Assumes a single HRF term in \code{event_model} and returns an
#' \code{hrfals_fit} object.
#'
#' @inheritParams fmrireg_hrf_cfals
#' @param baseline_model Optional baseline model to project out alongside
#'   `confound_obj`.
#' @param beta_penalty List with elements `l1`, `alpha`, and `warm_start`
#'   forwarded to `estimate_hrf_cfals` for sparse beta estimation.
#' @param design_control List of design matrix processing options. Set
#'   `standardize_predictors = TRUE` to z-score continuous predictors
#'   before estimation.
#' @return An object of class \code{hrfals_fit}.
#' @export
hrfals <- function(fmri_data_obj,
                   event_model,
                   hrf_basis,
                   confound_obj = NULL,
                   baseline_model = NULL,
                   beta_penalty = list(l1 = 0, alpha = 1, warm_start = TRUE),
                   design_control = list(standardize_predictors = TRUE,
                                         cache_design_blocks = TRUE),
                   lam_beta = 10,
                   lam_h = 1,
                   lambda_s = 0,
                   laplacian_obj = NULL,
                   h_solver = c("direct", "cg", "auto"),
                   cg_max_iter = 100,
                   cg_tol = 1e-4,
                   R_mat = NULL,
                   fullXtX = FALSE,
                   max_alt = 1,
                   ...) {
  target_term <- names(event_model$terms)[1]
  R_mat_param <- if (is.null(R_mat)) "identity" else R_mat
  estimate_hrf_cfals(fmri_data_obj,
                     fmrireg_event_model = event_model,
                     target_event_term_name = target_term,
                     hrf_basis_for_cfals = hrf_basis,
                     confound_obj = confound_obj,
                     baseline_model = baseline_model,
                     method = "cf_als",
                     lambda_b = lam_beta,
                     lambda_h = lam_h,
                     lambda_s = lambda_s,
                     laplacian_obj = laplacian_obj,
                     h_solver = h_solver,
                     cg_max_iter = cg_max_iter,
                     cg_tol = cg_tol,
                     R_mat = R_mat_param,
                     fullXtX = fullXtX,
                     max_alt = max_alt,
                     beta_penalty = beta_penalty,
                     design_control = design_control,
                     ...)
}

#' @keywords internal
fmrireg_cfals <- function(...) {
  .Deprecated("hrfals")
  hrfals(...)
}

#' @keywords internal
fmrireg_hrf_cfals <- function(...) {
  .Deprecated("hrfals")
  hrfals(...)
}
