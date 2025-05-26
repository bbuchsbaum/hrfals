#' LS+SVD Rank-1 HRF Estimation Engine
#'
#' Internal helper implementing the LS+SVD algorithm described in
#' `data-raw/LSS+SVD_proposal.md`. The design matrices should already
#' be projected to the confound null space.
#'
#' @param X_list_proj list of k design matrices (n x d each)
#' @param Y_proj      numeric matrix of projected BOLD data (n x v)
#' @param lambda_init ridge penalty for initial GLM solve
#' @param Phi_recon_matrix Reconstruction matrix mapping coefficients to HRF
#'   shape (p x d)
#' @param h_ref_shape_canonical Canonical reference HRF shape of length p for
#'   sign alignment
#' @param svd_backend currently ignored, placeholder for future backends
#' @param epsilon_svd tolerance for singular value screening
#' @param epsilon_scale tolerance for scale in identifiability step
#' @param R_mat optional penalty matrix for the initial GLM solve
#' @return list with matrices `h` (d x v), `beta` (k x v) and
#'         `Gamma_hat` (d*k x v)
#' @keywords internal
#' @noRd
ls_svd_engine <- function(X_list_proj, Y_proj, lambda_init = 1,
                          Phi_recon_matrix,
                          h_ref_shape_canonical,
                          svd_backend = c("base_R"),
                          epsilon_svd = 1e-8,
                          epsilon_scale = 1e-8,
                          R_mat = NULL) {
  svd_backend <- match.arg(svd_backend)
  stopifnot(is.list(X_list_proj), length(X_list_proj) >= 1)
  n <- nrow(Y_proj)
  v <- ncol(Y_proj)
  d <- ncol(X_list_proj[[1]])
  k <- length(X_list_proj)
  for (X in X_list_proj) {
    if (nrow(X) != n) stop("Design matrices must have same rows as Y_proj")
    if (ncol(X) != d) stop("All design matrices must have the same column count")
  }

  if (!is.matrix(Phi_recon_matrix) || ncol(Phi_recon_matrix) != d)
    stop("`Phi_recon_matrix` must be a p x d matrix")
  if (length(h_ref_shape_canonical) != nrow(Phi_recon_matrix))
    stop("`h_ref_shape_canonical` must have length nrow(Phi_recon_matrix)")

  if (!is.null(R_mat)) {
    if (!is.matrix(R_mat) || nrow(R_mat) != d || ncol(R_mat) != d) {
      stop(paste("R_mat must be a", d, "x", d, "matrix"))
    }
  }

  cholSolve <- function(M, B, eps = max(epsilon_svd, epsilon_scale)) {
    L <- tryCatch(chol(M),
                  error = function(e) chol(M + eps * diag(nrow(M))))
    backsolve(L, forwardsolve(t(L), B))
  }

  Xbig <- do.call(cbind, X_list_proj)
  XtX  <- crossprod(Xbig)
  Xty  <- crossprod(Xbig, Y_proj)
  penalty_mat <- if (is.null(R_mat)) diag(d) else R_mat
  R_big <- kronecker(diag(k), penalty_mat)
  XtX_ridge <- XtX + lambda_init * R_big
  Gamma_hat <- cholSolve(XtX_ridge, Xty)

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

  H_shapes <- Phi_recon_matrix %*% H_out
  scl <- apply(abs(H_shapes), 2, max)
  flip <- rep(1.0, v)
  align <- colSums(H_shapes * h_ref_shape_canonical)
  flip[align < 0 & scl > epsilon_scale] <- -1.0
  eff_scl <- pmax(scl, epsilon_scale)
  H_final <- sweep(H_out, 2, flip / eff_scl, "*")
  B_final <- sweep(B_out, 2, flip * eff_scl, "*")
  zero_idx <- scl <= epsilon_scale
  if (any(zero_idx)) {
    H_final[, zero_idx] <- 0
    B_final[, zero_idx] <- 0
  }

  dimnames(H_final) <- list(NULL, colnames(Y_proj))
  dimnames(B_final) <- list(names(X_list_proj), colnames(Y_proj))

  list(h = H_final, beta = B_final, Gamma_hat = Gamma_hat)
}

