#' LS+SVD+1ALS Rank-1 HRF Estimation Engine
#'
#' Internal helper implementing the LS+SVD+1ALS algorithm described in
#' `data-raw/LSS+SVD_proposal.md`. Inputs should be projected to the
#' confound null space.
#'
#' @param X_list_proj list of k design matrices (n x d each)
#' @param Y_proj numeric matrix of projected BOLD data (n x v)
#' @param lambda_init ridge penalty for initial GLM solve
#' @param lambda_b ridge penalty for \eqn{\beta}-update
#' @param lambda_h ridge penalty for \eqn{h}-update
#' @param lambda_joint joint penalty for \eqn{\beta} and \eqn{h} updates
#' @param fullXtX_flag Logical. If `TRUE`, the h-update step uses the full
#'   Gramian including cross-condition terms
#'   \eqn{(\sum_l \beta_l X_l)^\top(\sum_m \beta_m X_m)}. If `FALSE` (default)
#'   the Gramian omits cross-condition terms,
#'   \eqn{\sum_l \beta_l^2 X_l^\top X_l}. A single shared HRF coefficient
#'   vector is still estimated per voxel.
#' @param Phi_recon_matrix Reconstruction matrix mapping coefficients to HRF
#'   shape (p x d)
#' @param h_ref_shape_canonical Canonical reference HRF shape of length p for
#'   sign alignment
#' @param svd_backend backend for SVD in the initialization step
#' @param epsilon_svd tolerance for singular value screening
#' @param epsilon_scale tolerance for scale in identifiability step
#' @param R_mat optional penalty matrix for the h-update
#' @return list with matrices `h` (d x v), `beta` (k x v), `Gamma_hat` (d*k x v)
#'         and the initial estimates `h_ls_svd`, `beta_ls_svd`
#' @keywords internal
#' @noRd
ls_svd_1als_engine <- function(X_list_proj, Y_proj,
                               lambda_init = 1,
                               lambda_b = 10,
                               lambda_h = 1,
                               lambda_joint = 0,
                               fullXtX_flag = FALSE,
                               Phi_recon_matrix,
                               h_ref_shape_canonical,
                               svd_backend = c("base_R"),
                               epsilon_svd = 1e-8,
                               epsilon_scale = 1e-8,
                               R_mat = NULL) {

  if (lambda_init < 0 || lambda_b < 0 || lambda_h < 0 || lambda_joint < 0)
    stop("Lambdas must be non-negative")
  
  # Validate inputs and extract dimensions
  dims <- validate_hrf_engine_inputs(X_list_proj, Y_proj, Phi_recon_matrix, h_ref_shape_canonical)
  d <- dims$d

  if (!is.null(R_mat)) {
    if (!is.matrix(R_mat) || nrow(R_mat) != d || ncol(R_mat) != d) {
      stop(paste("R_mat must be a", d, "x", d, "matrix"))
    }
  }

  init <- ls_svd_engine(X_list_proj, Y_proj,
                        lambda_init = lambda_init,
                        Phi_recon_matrix = Phi_recon_matrix,
                        h_ref_shape_canonical = h_ref_shape_canonical,
                        svd_backend = svd_backend,
                        epsilon_svd = epsilon_svd,
                        epsilon_scale = epsilon_scale)

  h_current <- init$h
  b_current <- init$beta
  k <- dims$k
  v <- dims$v

  XtX_list <- lapply(X_list_proj, crossprod)
  XtY_list <- lapply(X_list_proj, function(X) crossprod(X, Y_proj))

  XtX_full_list <- NULL
  if (fullXtX_flag) {
    XtX_full_list <- matrix(vector("list", k * k), k, k)
    for (l in seq_len(k)) {
      for (m in seq_len(k)) {
        XtX_full_list[[l, m]] <- crossprod(X_list_proj[[l]], X_list_proj[[m]])
      }
    }
  }

  H_als <- matrix(0.0, d, v)
  B_als <- matrix(0.0, k, v)

  for (vx in seq_len(v)) {
    h_vx <- h_current[, vx]
    DhTy_vx <- vapply(seq_len(k), function(c)
      crossprod(h_vx, XtY_list[[c]][, vx]), numeric(1))
    if (fullXtX_flag) {
      G_vx <- matrix(0.0, k, k)
      for (l in seq_len(k)) {
        for (m in seq_len(k)) {
          G_vx[l, m] <- crossprod(h_vx, XtX_full_list[[l, m]] %*% h_vx)
        }
      }
    } else {
      diag_vals <- vapply(seq_len(k), function(c)
        crossprod(h_vx, XtX_list[[c]] %*% h_vx), numeric(1))
      G_vx <- diag(diag_vals, k)
    }
      # Apply joint ridge penalty to beta update
      G_vx <- G_vx + lambda_joint * diag(k)
      B_als[, vx] <- cholSolve(G_vx + lambda_b * diag(k), DhTy_vx,
                               eps = max(epsilon_svd, epsilon_scale))
  }

  b_current <- B_als

  for (vx in seq_len(v)) {
    b_vx <- b_current[, vx]
    penalty_mat <- if (is.null(R_mat)) diag(d) else R_mat
    lhs <- lambda_h * penalty_mat + lambda_joint * diag(d)
    rhs <- numeric(d)
    for (l in seq_len(k)) {
      rhs <- rhs + b_vx[l] * XtY_list[[l]][, vx]
      if (fullXtX_flag) {
        for (m in seq_len(k)) {
          lhs <- lhs + b_vx[l] * b_vx[m] * XtX_full_list[[l, m]]
        }
      } else {
        lhs <- lhs + b_vx[l]^2 * XtX_list[[l]]
      }
    }
      H_als[, vx] <- cholSolve(lhs, rhs,
                               eps = max(epsilon_svd, epsilon_scale))
  }

  # Normalize and align HRF shapes
  result <- normalize_and_align_hrf(H_als, B_als, Phi_recon_matrix, 
                                   h_ref_shape_canonical, epsilon_scale,
                                   Y_proj, X_list_proj)

  list(h = result$h, beta = result$beta,
       h_ls_svd = init$h, beta_ls_svd = init$beta, Gamma_hat = init$Gamma_hat)
}

