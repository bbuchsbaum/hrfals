#' Fit Rank-1 HRF Using CF-ALS
#'
#' Convenience wrapper for the original CF-ALS implementation. This
#' function simply calls [estimate_hrf_cfals()] with CF-ALS defaults and
#' retains the historical argument names for backwards compatibility.
#'
#' @param fmri_data_obj An `fmri_dataset` or numeric matrix of BOLD data
#'   (time points x voxels). If a dataset, sampling information is
#'   taken from the object.
#' @param event_model An `event_model` describing the stimuli to use
#'   for HRF estimation.
#' @param hrf_basis An `HRF` basis object from the `fmrihrf` package
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
#' This function implements the Confound-Adjusted Alternating Least Squares (CF-ALS)
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
#'   \item SVD Initialization: Robust starting values for `h` and `beta` are obtained using a regularized least squares solution followed by Singular Value Decomposition (SVD) of the initial coefficient estimates.
#'   \item CF-ALS Alternation: The algorithm alternates between updating `beta` (amplitudes) holding `h` (HRF coefficients) fixed, and updating `h` holding `beta` fixed. This process is repeated for `max_alt` iterations.
#'     \itemize{
#'       \item β-update: Solves (G_mat + λ_β I)β_v = D(h)^T y_v for each voxel v, where G_mat = D(h)^T D(h).
#'       \item h-update: Solves (LHS + λ_h I)h = RHS, where LHS and RHS depend on the data and current beta estimates. The `fullXtX` parameter influences the construction of LHS.
#'     }
#'   \item Identifiability: The estimated HRF coefficients `h` are typically normalized (e.g., ||Φh||_∞ = 1, where Φ is the basis reconstruction matrix) and their sign aligned with a canonical HRF shape to ensure consistent scaling and polarity across voxels/conditions.
#' }
#'
#' The function is designed to be compatible with any HRF basis defined in the
#' `fmrihrf` package (provided `nbasis(hrf_basis) > 1`). The `R_mat`
#' parameter allows for basis-specific penalty matrices for HRF coefficient
#' regularization, although the default is an identity matrix (standard ridge penalty).
#'
#' R\eqn{^2} (coefficient of determination) is computed on the data *after*
#' any confound projection has been applied.
#'
#' @seealso [fmrireg_cfals()] for the deprecated alias, [estimate_hrf_cfals()] for the lower-level interface.
#' @references (If applicable, add references to papers describing the CF-ALS method or its implementation).
#'
#' @examples
#' # Simulate data
#' sframe <- fmridesign::sampling_frame(blocklens = 40, TR = 1)
#' ev_df <- data.frame(onset = c(5, 15, 25), block = 1, cond = "A")
#' emod <- fmridesign::event_model(onset ~ fmridesign::hrf(cond, basis = fmrihrf::HRF_SPMG3),
#'                              data = ev_df, block = ~ block,
#'                              sampling_frame = sframe)
#' Y_matrix <- matrix(rnorm(40 * 5), 40, 5) # 40 timepoints, 5 voxels
#'
#' # Fit using hrfals CF-ALS defaults
#' cfals_fit <- hrfals(
#'   fmri_data_obj = Y_matrix,
#'   event_model = emod,
#'   hrf_basis = fmrihrf::HRF_SPMG3, # Using SPMG3 basis (3 basis functions)
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
#' bspline_basis <- fmrihrf::HRF_BSPLINE
#' emod_bspline <- fmridesign::event_model(onset ~ fmridesign::hrf(cond, basis = bspline_basis),
#'                                      data = ev_df, block = ~ block,
#'                                      sampling_frame = sframe)
#' cfals_fit_bspline <- hrfals(
#'   fmri_data_obj = Y_matrix,
#'   event_model = emod_bspline,
#'   hrf_basis = bspline_basis,
#'   max_alt = 2
#' )
#' print(cfals_fit_bspline)
#' @param fmri_data_obj An \code{fmri_dataset} object or numeric matrix of 
#'   BOLD data (time points x voxels).
#' @param event_model An \code{event_model} object describing the experimental
#'   design and stimulus events.
#' @param hrf_basis An \code{HRF} basis object from the fmrihrf package
#'   defining the basis functions for HRF estimation.
#' @param confound_obj Optional matrix of confound regressors to project out
#'   (e.g., motion parameters, physiological noise).
#' @param baseline_model Optional baseline model to project out alongside
#'   \code{confound_obj}.
#' @param beta_penalty List with elements \code{l1}, \code{alpha}, and 
#'   \code{warm_start} for sparse beta estimation. Setting \code{l1 > 0} 
#'   triggers an Elastic Net beta update with mixing parameter \code{alpha}. 
#'   Typical \code{l1} values are between 0.01 and 0.1 for scaled predictors. 
#'   Warm-starting reuses betas from the previous iteration to speed convergence.
#' @param design_control List of design matrix processing options. Set
#'   \code{standardize_predictors = TRUE} to z-score continuous predictors
#'   before estimation. Betas are returned in the original units. When
#'   \code{cache_design_blocks = TRUE}, lagged design matrices are cached in
#'   memory when feasible to avoid recomputation.
#' @param lam_beta Ridge penalty for the beta (amplitude) update step.
#' @param lam_h Ridge penalty for the h (HRF shape) update step.
#' @param lambda_s Spatial regularization strength for smoothing HRF estimates
#'   across neighboring voxels. Requires \code{laplacian_obj} when > 0.
#' @param laplacian_obj Optional Laplacian matrix object from 
#'   \code{build_voxel_laplacian()} for spatial regularization.
#' @param h_solver Solver for the h-update: "direct" (default), "cg" 
#'   (conjugate gradient), or "auto".
#' @param cg_max_iter Maximum iterations for conjugate gradient solver.
#' @param cg_tol Convergence tolerance for conjugate gradient solver.
#' @param R_mat Optional penalty matrix for HRF coefficients. If NULL,
#'   uses identity matrix (standard ridge penalty).
#' @param fullXtX Logical. If TRUE, uses full cross-condition Gramian matrix
#'   in h-update. If FALSE (default), uses diagonal approximation.
#' @param max_alt Maximum number of alternating least squares iterations.
#' @param ... Additional arguments passed to \code{estimate_hrf_cfals}.
#' 
#' @return An object of class \code{hrfals_fit} containing:
#'   \item{h_coeffs}{Matrix of HRF basis coefficients (nbasis x voxels)}
#'   \item{beta_amps}{Matrix of condition amplitudes (conditions x voxels)}
#'   \item{recon_hrf}{Reconstructed HRF time courses}
#'   \item{gof}{Goodness-of-fit (R-squared) per voxel}
#'   \item{residuals}{Model residuals}
#'   \item{method}{Estimation method used}
#'   \item{lambdas}{List of regularization parameters used}
#'   
#' @examples
#' \donttest{
#' library(fmridesign)
#' # Create sampling frame and event model
#' sf <- sampling_frame(blocklens = 60, TR = 1)
#' events <- data.frame(
#'   onset = c(5, 15, 30, 45),
#'   condition = factor(c("A", "A", "B", "B")),
#'   block = 1
#' )
#' emod <- event_model(onset ~ hrf(condition), data = events,
#'                     block = ~ block, sampling_frame = sf)
#' 
#' # Simulate BOLD data
#' Y <- matrix(rnorm(60 * 10), 60, 10)  # 60 timepoints, 10 voxels
#' 
#' # Fit HRF using CF-ALS
#' fit <- hrfals(Y, emod, fmrihrf::HRF_SPMG3)
#' print(fit)
#' }
#' 
#' @seealso \code{\link{hrfals_sparse}} for sparse estimation,
#'   \code{\link{estimate_hrf_cfals}} for lower-level interface,
#'   \code{\link{create_cfals_design}} for design matrix creation
#'   
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
                     fmridesign_event_model = event_model,
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
  .Deprecated("hrfals",
              msg = "`fmrireg_cfals()` is deprecated; use `hrfals()` instead.")
  hrfals(...)
}

#' @keywords internal
fmrireg_hrf_cfals <- function(...) {
  .Deprecated("hrfals",
              msg = "`fmrireg_hrf_cfals()` is deprecated; use `hrfals()` instead.")
  hrfals(...)
}

#' Convenient wrapper enabling sparse beta estimation
#'
#' `hrfals_sparse()` simply calls [hrfals()] with sparse beta defaults
#' (`beta_penalty$l1 > 0`) and standard design control. It provides a
#' shorthand for common sparse modeling use cases.
#'
#' @inheritParams hrfals
#' @param ... Arguments passed on to [hrfals()].
#' @export
hrfals_sparse <- function(...,
                          beta_penalty = list(l1 = 0.05,
                                               alpha = 1,
                                               warm_start = TRUE),
                          design_control = list(standardize_predictors = TRUE,
                                                cache_design_blocks = TRUE)) {
  hrfals(...,
         beta_penalty = beta_penalty,
         design_control = design_control)
}
