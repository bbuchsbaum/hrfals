#' LS+SVD Rank-1 HRF Estimation Engine
#'
#' Internal helper implementing the LS+SVD algorithm described in
#' `data-raw/LSS+SVD_proposal.md`. The design matrices should already
#' be projected to the confound null space.
#'
#' @param X_list_proj list of k design matrices (n x d each)
#' @param Y_proj      numeric matrix of projected BOLD data (n x v)
#' @param lambda_init ridge penalty for initial GLM solve
#' @param Phi_recon_matrix Canonical reference HRF shape of length p for
#'   sign alignment
#' @param h_ref_shape_canonical Canonical reference HRF shape of length p for
#'   sign alignment
#' @param svd_backend currently ignored, placeholder for future backends
#' @param epsilon_svd tolerance for singular value screening
#' @param epsilon_scale tolerance for scale in identifiability step
#' @return list with matrices `h` (d x v), `beta` (k x v) and
#'         `Gamma_hat` (d*k x v)
#' @keywords internal
#' @noRd
ls_svd_engine <- function(X_list_proj, Y_proj, lambda_init = 1,
                          Phi_recon_matrix,
                          h_ref_shape_canonical,
                          svd_backend = c("base_R"),
                          epsilon_svd = 1e-8,
                          epsilon_scale = 1e-8) {
  svd_backend <- match.arg(svd_backend)
  
  # Validate inputs and extract dimensions
  dims <- validate_hrf_engine_inputs(X_list_proj, Y_proj, Phi_recon_matrix, h_ref_shape_canonical)
  n <- dims$n
  v <- dims$v
  d <- dims$d
  k <- dims$k


  Xbig <- do.call(cbind, X_list_proj)
  XtX  <- crossprod(Xbig)
  Xty  <- crossprod(Xbig, Y_proj)
  XtX_ridge <- XtX + lambda_init * diag(d * k)
  Gamma_hat <- cholSolve(XtX_ridge, Xty, eps = max(epsilon_svd, epsilon_scale))

  H_out <- matrix(0.0, d, v)
  B_out <- matrix(0.0, k, v)

  for (vx in seq_len(v)) {
    G_vx <- matrix(Gamma_hat[, vx], nrow = d, ncol = k)
    sv <- svd(G_vx, nu = 1, nv = 1)
    if (length(sv$d) && sv$d[1] > epsilon_svd) {
      s1 <- sqrt(sv$d[1])
      H_out[, vx] <- sv$u[, 1] * s1
      B_out[, vx] <- sv$v[, 1] * s1
    }
  }

  # Normalize and align HRF shapes
  result <- normalize_and_align_hrf(H_out, B_out, Phi_recon_matrix, 
                                   h_ref_shape_canonical, epsilon_scale,
                                   Y_proj, X_list_proj)

  list(h = result$h, beta = result$beta, Gamma_hat = Gamma_hat)
}

