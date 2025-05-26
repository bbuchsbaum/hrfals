#' Confound-Free ALS HRF Estimation Engine
#'
#' Internal helper implementing the CF-ALS algorithm described in
#' `data-raw/CFALS_proposal.md`. Inputs should already be projected
#' to the confound null space.
#'
#' @param X_list_proj list of k design matrices (n x d each)
#' @param Y_proj numeric matrix of projected BOLD data (n x v)
#' @param lambda_b ridge penalty for beta-update
#' @param lambda_h ridge penalty for h-update
#' @param R_mat_eff effective d x d penalty matrix for h-update. If NULL a
#'   diagonal matrix is used.
#' @param fullXtX_flag logical; if TRUE use cross-condition terms in h-update
#' @param precompute_xty_flag logical; if TRUE precompute `XtY_list` otherwise
#'   compute per voxel on-the-fly
#' @param Phi_recon_matrix Reconstruction matrix mapping coefficients to HRF
#'   shape (p x d)
#' @param h_ref_shape_canonical Canonical reference HRF shape of length p for
#'   sign alignment
#' @param max_alt number of alternating updates after initialization
#' @param epsilon_svd tolerance for singular value screening
#' @param epsilon_scale tolerance for scale in identifiability step
#' @return list with matrices `h` (d x v) and `beta` (k x v). The
#'   matrix `h` has an attribute `"iterations"` recording the number
#'   of alternating updates performed.
#' @keywords internal
#' @noRd
cf_als_engine <- function(X_list_proj, Y_proj,
                          lambda_b = 10,
                          lambda_h = 1,
                          R_mat_eff = NULL,
                          fullXtX_flag = FALSE,
                          precompute_xty_flag = TRUE,

                          Phi_recon_matrix,
                          h_ref_shape_canonical,
                          max_alt = 1,
                          epsilon_svd = 1e-8,
                          epsilon_scale = 1e-8,
                          Phi_recon_matrix = NULL,
                          h_ref_shape_canonical = NULL) {

  stopifnot(is.list(X_list_proj), length(X_list_proj) >= 1)
  n <- nrow(Y_proj)
  v <- ncol(Y_proj)
  d <- ncol(X_list_proj[[1]])
  k <- length(X_list_proj)
  if (!is.matrix(Phi_recon_matrix) || ncol(Phi_recon_matrix) != d)
    stop("`Phi_recon_matrix` must be a p x d matrix")
  if (length(h_ref_shape_canonical) != nrow(Phi_recon_matrix))
    stop("`h_ref_shape_canonical` must have length nrow(Phi_recon_matrix)")
  for (X in X_list_proj) {
    if (nrow(X) != n) stop("Design matrices must have same rows as Y_proj")
    if (ncol(X) != d) stop("All design matrices must have the same column count")
  }
  if (lambda_b < 0 || lambda_h < 0) {
    stop("lambda_b and lambda_h must be non-negative")
  }

  if (!is.null(R_mat_eff)) {
    if (!is.matrix(R_mat_eff) || nrow(R_mat_eff) != d || ncol(R_mat_eff) != d) {
      stop(paste("R_mat_eff must be a d x d matrix, where d is", d))
    }
  }

  cholSolve <- function(M, b, eps = max(epsilon_svd, epsilon_scale)) {
    ok <- TRUE
    L  <- tryCatch(chol(M), error = function(e) { ok <<- FALSE ; NULL })
    if (!ok || (is.null(L) || min(diag(L)) < eps)) {
      L <- chol(M + eps * diag(nrow(M)))
    }
    backsolve(L, forwardsolve(t(L), b))
  }

  init <- ls_svd_engine(X_list_proj, Y_proj,
                        lambda_init = 0,
                        Phi_recon_matrix = Phi_recon_matrix,
                        h_ref_shape_canonical = h_ref_shape_canonical,
                        epsilon_svd = epsilon_svd,
                        epsilon_scale = epsilon_scale)
  h_current <- init$h
  b_current <- init$beta

  XtX_list <- lapply(X_list_proj, crossprod)


  size_est <- as.numeric(k) * d * v * 8
  if (precompute_xty_flag && size_est > 2e9) {
    message("Estimated size of XtY_list (", size_est,
            " bytes) is large; consider `precompute_xty_flag = FALSE`")
  }
  XtY_list <- if (precompute_xty_flag) {
    lapply(X_list_proj, function(X) crossprod(X, Y_proj))
  } else {
    NULL

  }

  if (fullXtX_flag) {
    XtX_full_list <- matrix(vector("list", k * k), k, k)
    for (l in seq_len(k)) {
      for (m in seq_len(k)) {
        XtX_full_list[[l, m]] <- crossprod(X_list_proj[[l]], X_list_proj[[m]])
      }
    }
  } else {
    XtX_full_list <- NULL
  }

  iter_final <- 0
  for (iter in seq_len(max_alt)) {
    b_prev <- b_current
    h_prev <- h_current
    for (vx in seq_len(v)) {
      h_vx <- h_current[, vx]

      if (!isTRUE(precompute_xty_flag)) {
        XtY_cache <- vector("list", k)
        for (l in seq_len(k)) {
          XtY_cache[[l]] <- crossprod(X_list_proj[[l]], Y_proj[, vx])
        }
      }

      DhTy_vx <- vapply(seq_len(k), function(c) {
        XtY_c_vx <- if (isTRUE(precompute_xty_flag)) {
          XtY_list[[c]][, vx]
        } else {
          XtY_cache[[c]]
        }
        crossprod(h_vx, XtY_c_vx)

      }, numeric(1))
      G_vx <- matrix(0.0, k, k)
      for (l in seq_len(k)) {
        if (fullXtX_flag) {
          for (m in seq_len(k)) {
            term <- XtX_full_list[[l, m]]
            G_vx[l, m] <- crossprod(h_vx, term %*% h_vx)
          }
        } else {
          G_vx[l, l] <- crossprod(h_vx, XtX_list[[l]] %*% h_vx)
        }
      }
      b_current[, vx] <- cholSolve(G_vx + lambda_b * diag(k), DhTy_vx)
    }

    h_penalty_matrix <- if (is.null(R_mat_eff)) {
      diag(d)
    } else {
      Matrix::forceSymmetric(R_mat_eff)
    }

    for (vx in seq_len(v)) {
      b_vx <- b_current[, vx]
      lhs <- lambda_h * h_penalty_matrix
      rhs <- numeric(d)
      for (l in seq_len(k)) {

        XtY_l_vx <- if (isTRUE(precompute_xty_flag)) {
          XtY_list[[l]][, vx]
        } else {
          XtY_cache[[l]]
        }
        rhs <- rhs + b_vx[l] * XtY_l_vx
        if (fullXtX_flag) {
          for (m in seq_len(k)) {
            lhs <- lhs + b_vx[l] * b_vx[m] * XtX_full_list[[l, m]]
          }
        } else {
          lhs <- lhs + b_vx[l]^2 * XtX_list[[l]]
        }
      }
    h_current[, vx] <- cholSolve(lhs, rhs)
    }

    iter_final <- iter
    if (max(abs(b_current - b_prev)) < 1e-6 &&
        max(abs(h_current - h_prev)) < 1e-6) {
      break
    }
  }

  H_shapes_iter <- Phi_recon_matrix %*% h_current
  scl <- apply(abs(H_shapes_iter), 2, max)
  flip <- rep(1.0, v)
  align_scores <- colSums(H_shapes_iter * h_ref_shape_canonical)
  flip[align_scores < 0 & scl > epsilon_scale] <- -1.0
  eff_scl <- pmax(scl, epsilon_scale)
  h_final <- sweep(h_current, 2, flip / eff_scl, "*")
  b_final <- sweep(b_current, 2, flip * eff_scl, "*")
  zero_idx <- scl <= epsilon_scale
  if (any(zero_idx)) {
    h_final[, zero_idx] <- 0
    b_final[, zero_idx] <- 0
  }

  attr(h_final, "iterations") <- iter_final
  list(h = h_final, beta = b_final)
}

