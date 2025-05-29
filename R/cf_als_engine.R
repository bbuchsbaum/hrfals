#' Compute LHS Blocks for the h-update
#'
#' Internal helper used by `cf_als_engine` to construct the list of
#' `d x d` left-hand-side matrices for each voxel's h-update step.
#'
#' @param XtX_list list of `k` cross-product matrices `X_l^T X_l`.
#' @param XtX_full_list optional `k x k` list matrix of full cross-products
#'   `X_l^T X_m` when `fullXtX_flag` is `TRUE`.
#' @param b_current current beta estimates (k x v).
#' @param lambda_h ridge penalty for the h-update.
#' @param lambda_joint joint ridge penalty applied to both beta and h.
#' @param R_mat_eff effective penalty matrix for h-update or `NULL`.
#' @param fullXtX_flag logical; whether full Gramian is used.
#' @param d number of HRF basis functions.
#' @param v number of voxels.
#' @param k number of conditions.
#' @return list of length `v` with `d x d` matrices.
#' @keywords internal
#' @noRd
make_lhs_block_list <- function(XtX_list, XtX_full_list, b_current,
                                lambda_h, lambda_joint, R_mat_eff,
                                fullXtX_flag, d, v, k) {
  penalty_mat <- if (is.null(R_mat_eff)) diag(d) else R_mat_eff
  base_lhs <- lambda_h * penalty_mat + lambda_joint * diag(d)
  lhs_list <- vector("list", v)
  for (vx in seq_len(v)) {
    lhs_vx <- base_lhs
    b_vx <- b_current[, vx]
    if (fullXtX_flag) {
      for (l in seq_len(k)) {
        for (m in seq_len(k)) {
          lhs_vx <- lhs_vx + b_vx[l] * b_vx[m] * XtX_full_list[[l, m]]
        }
      }
    } else {
      for (l in seq_len(k)) {
        lhs_vx <- lhs_vx + b_vx[l]^2 * XtX_list[[l]]
      }
    }
    lhs_list[[vx]] <- lhs_vx
  }
  lhs_list
}

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
#' @param lambda_init ridge penalty for initial LS+SVD step
#' @param lambda_joint joint penalty for beta-h update
#' @param R_mat_eff effective d x d penalty matrix for h-update. If NULL a
#'   diagonal matrix is used.
#' @param fullXtX_flag Logical. If `TRUE`, the h-update step uses the full
#'   Gramian including cross-condition terms
#'   \eqn{(\sum_l \beta_l X_l)^\top(\sum_m \beta_m X_m)}. If `FALSE`
#'   (default) the Gramian omits cross-condition terms,
#'   \eqn{\sum_l \beta_l^2 X_l^\top X_l}. A single shared HRF coefficient
#'   vector is still estimated per voxel.
#' @param precompute_xty_flag logical; if TRUE precompute `XtY_list` otherwise
#'   compute per voxel on-the-fly
#' @param Phi_recon_matrix p x d matrix for sign alignment
#' @param h_ref_shape_canonical Canonical reference HRF shape of length p for
#'   sign alignment
#' @param max_alt number of alternating updates after initialization
#' @param epsilon_svd tolerance for singular value screening
#' @param epsilon_scale tolerance for scale in identifiability step
#' @param beta_penalty List with elements `l1`, `alpha`, and `warm_start`
#'   controlling sparse (Elastic Net) penalties for the beta-step.
#'   Set `l1 > 0` to enable sparse estimation.
#' @param design_control List of design matrix processing options. Set
#'   `standardize_predictors = TRUE` to z-score predictors before estimation
#'   and `cache_design_blocks = TRUE` to pre-compute lagged design blocks.
#' @return list with matrices `h` (d x v) and `beta` (k x v). The
#'   matrix `h` has an attribute `"iterations"` recording the number
#'   of alternating updates performed.
#' @keywords internal
#' @importFrom Matrix forceSymmetric
#' @noRd
cf_als_engine <- function(X_list_proj, Y_proj,
                          lambda_b = 10,
                          lambda_h = 1,
                          lambda_init = 1,
                          lambda_joint = 0,
                          R_mat_eff = NULL,
                          fullXtX_flag = FALSE,
                          precompute_xty_flag = TRUE,
                          lambda_s = 0,
                          laplacian_obj = NULL,
                          h_solver = c("direct", "cg", "auto"),
                          cg_max_iter = 100,
                          cg_tol = 1e-4,
                          Phi_recon_matrix,
                          h_ref_shape_canonical,
                          max_alt = 1,
                          epsilon_svd = 1e-8,
                          epsilon_scale = 1e-8,
                          beta_penalty = list(l1 = 0, alpha = 1, warm_start = TRUE),
                          design_control = list(standardize_predictors = TRUE,
                                                cache_design_blocks = TRUE)) {

  # Validate inputs and extract dimensions
  dims <- validate_hrf_engine_inputs(X_list_proj, Y_proj, Phi_recon_matrix, h_ref_shape_canonical)
  n <- dims$n
  v <- dims$v
  d <- dims$d
  k <- dims$k
  if (lambda_b < 0 || lambda_h < 0 || lambda_init < 0 || lambda_joint < 0) {
    stop("lambda_b, lambda_h, lambda_init, and lambda_joint must be non-negative")
  }

  h_solver <- match.arg(h_solver)
  if (lambda_s < 0) {
    stop("lambda_s must be non-negative")
  }

  if (lambda_s > 0) {
    if (is.null(laplacian_obj) || is.null(laplacian_obj$L) ||
        is.null(laplacian_obj$degree)) {
      stop("laplacian_obj with elements L and degree must be provided when lambda_s > 0")
    }
    L_mat <- laplacian_obj$L
    degree_vec <- laplacian_obj$degree
    if (length(degree_vec) != v) {
      stop("laplacian_obj$degree length mismatch with number of voxels")
    }
    if (h_solver == "auto") {
      current_solver <- if (d * v < 50000) "direct" else "cg"
    } else {
      current_solver <- h_solver
    }
  } else {
    L_mat <- NULL
    degree_vec <- NULL
    current_solver <- "direct"
  }

  if (!is.null(R_mat_eff)) {
    if (!is.matrix(R_mat_eff) || nrow(R_mat_eff) != d || ncol(R_mat_eff) != d) {
      stop(paste("R_mat_eff must be a d x d matrix, where d is", d))
    }
    # Force penalty matrix to be symmetric for numerical stability
    R_mat_eff <- as.matrix(Matrix::forceSymmetric(R_mat_eff))
  }


  init <- ls_svd_engine(X_list_proj, Y_proj,
                        lambda_init = lambda_init,
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
  obj_trace <- numeric(max_alt)
  perform_warm_start <- isTRUE(beta_penalty$warm_start) && beta_penalty$l1 > 0
  for (iter in seq_len(max_alt)) {
    b_prev <- b_current
    h_prev <- h_current

    if (beta_penalty$l1 > 0) {
      Xh_list <- lapply(X_list_proj, function(X) X %*% h_current)
      if (perform_warm_start && iter == 1) {
        warm_beta <- matrix(0.0, k, v)
        for (vx in seq_len(v)) {
          if (!isTRUE(precompute_xty_flag)) {
            XtY_cache <- lapply(X_list_proj, function(X) crossprod(X, Y_proj[, vx]))
          }

          h_vx <- h_current[, vx]
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
          G_vx <- G_vx + lambda_joint * diag(k)
          warm_beta[, vx] <- cholSolve(G_vx + lambda_b * diag(k), DhTy_vx,
                                       eps = max(epsilon_svd, epsilon_scale))
        }
        b_current <- warm_beta
      }
    }

    for (vx in seq_len(v)) {
      if (!isTRUE(precompute_xty_flag)) {
        XtY_cache <- lapply(X_list_proj, function(X) crossprod(X, Y_proj[, vx]))
      }

      h_vx <- h_current[, vx]

      if (beta_penalty$l1 > 0) {
        Xh_mat <- vapply(seq_len(k), function(c) Xh_list[[c]][, vx], numeric(n))
        y_vx <- Y_proj[, vx]
        fit <- glmnet::glmnet(x = Xh_mat, y = y_vx,
                              alpha = beta_penalty$alpha,
                              lambda = beta_penalty$l1,
                              standardize = FALSE,
                              intercept = FALSE)
        b_current[, vx] <- as.numeric(glmnet::coef(fit, s = beta_penalty$l1))[-1]
      } else {
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
        G_vx <- G_vx + lambda_joint * diag(k)
        b_current[, vx] <- cholSolve(G_vx + lambda_b * diag(k), DhTy_vx,
                                     eps = max(epsilon_svd, epsilon_scale))
      }
    }

    lhs_block_list <- make_lhs_block_list(
      XtX_list, XtX_full_list, b_current,
      lambda_h, lambda_joint, R_mat_eff,
      fullXtX_flag, d, v, k
    )

    if (lambda_s > 0 && current_solver == "cg") {
      RHS_mat <- matrix(0.0, d, v)
      for (vx in seq_len(v)) {
        if (!isTRUE(precompute_xty_flag)) {
          XtY_cache <- lapply(X_list_proj, function(X) crossprod(X, Y_proj[, vx]))
        }
        b_vx <- b_current[, vx]
        rhs <- numeric(d)
        for (l in seq_len(k)) {
          XtY_l_vx <- if (isTRUE(precompute_xty_flag)) {
            XtY_list[[l]][, vx]
          } else {
            XtY_cache[[l]]
          }
          rhs <- rhs + b_vx[l] * XtY_l_vx
        }
        RHS_mat[, vx] <- rhs
      }

      A_H_spatial <- construct_A_H_sparse(lhs_block_list, lambda_s, L_mat, d, v)
      P_sparse <- construct_preconditioner(lhs_block_list, lambda_s, degree_vec, d, v)
      cg_sol <- Rlinsolve::lsolve.cg(A = A_H_spatial,
                                     B = as.vector(RHS_mat),
                                     preconditioner = P_sparse,
                                     xinit = as.vector(h_current),
                                     reltol = cg_tol,
                                     maxiter = cg_max_iter,
                                     adjsym = TRUE, verbose = FALSE)
      if (cg_sol$iter == cg_max_iter &&
          cg_sol$errors[length(cg_sol$errors)] > cg_tol) {
        warning("CG solver did not converge within max_iter for h-update.")
      }
      h_current <- matrix(cg_sol$x, d, v)
      for (vx in seq_len(v)) {
        s <- max(abs(h_current[, vx]), epsilon_scale)
        h_current[, vx] <- h_current[, vx] / s
        b_current[, vx] <- b_current[, vx] * s
      }
    } else if (lambda_s > 0 && current_solver == "direct") {
      RHS_mat <- matrix(0.0, d, v)
      for (vx in seq_len(v)) {
        if (!isTRUE(precompute_xty_flag)) {
          XtY_cache <- lapply(X_list_proj, function(X) crossprod(X, Y_proj[, vx]))
        }
        b_vx <- b_current[, vx]
        rhs <- numeric(d)
        for (l in seq_len(k)) {
          XtY_l_vx <- if (isTRUE(precompute_xty_flag)) {
            XtY_list[[l]][, vx]
          } else {
            XtY_cache[[l]]
          }
          rhs <- rhs + b_vx[l] * XtY_l_vx
        }
        RHS_mat[, vx] <- rhs
      }

      A_H_spatial <- construct_A_H_sparse(lhs_block_list, lambda_s, L_mat, d, v)
      h_vec <- Matrix::solve(A_H_spatial, as.vector(RHS_mat))
      h_current <- matrix(as.numeric(h_vec), d, v)
      for (vx in seq_len(v)) {
        s <- max(abs(h_current[, vx]), epsilon_scale)
        h_current[, vx] <- h_current[, vx] / s
        b_current[, vx] <- b_current[, vx] * s
      }
    } else {
      for (vx in seq_len(v)) {
        if (!isTRUE(precompute_xty_flag)) {
          XtY_cache <- lapply(X_list_proj, function(X) crossprod(X, Y_proj[, vx]))
        }
        b_vx <- b_current[, vx]
        lhs <- lhs_block_list[[vx]]
        rhs <- numeric(d)
        for (l in seq_len(k)) {
          XtY_l_vx <- if (isTRUE(precompute_xty_flag)) {
            XtY_list[[l]][, vx]
          } else {
            XtY_cache[[l]]
          }
          rhs <- rhs + b_vx[l] * XtY_l_vx
        }
        # Block solve for entire h vector for this voxel
        h_current[, vx] <- cholSolve(lhs, rhs,
                                     eps = max(epsilon_svd, epsilon_scale))
        s <- max(abs(h_current[, vx]), epsilon_scale)
        h_current[, vx] <- h_current[, vx] / s
        b_current[, vx] <- b_current[, vx] * s
      }
    }


    # compute penalised SSE objective for monitoring
    pred_iter <- matrix(0, n, v)
    for (c in seq_len(k)) {
      pred_iter <- pred_iter + (X_list_proj[[c]] %*% h_current) *
        matrix(rep(b_current[c, ], each = n), n, v)
    }
    res_iter <- Y_proj - pred_iter
    SSE_iter <- sum(res_iter^2)
    beta_pen <- lambda_b * sum(b_current^2)
    R_eff <- if (is.null(R_mat_eff)) diag(d) else R_mat_eff
    h_pen <- lambda_h * sum(colSums(h_current * (R_eff %*% h_current)))
    # Include joint ridge penalty in objective
    joint_pen <- lambda_joint * (sum(b_current^2) + sum(h_current^2))
    obj_trace[iter] <- SSE_iter + beta_pen + h_pen + joint_pen

    iter_final <- iter
    
    # Check convergence based on parameter changes
    # After scale normalization, this should work better
    if (iter > 1) {
      b_change <- max(abs(b_current - b_prev))
      h_change <- max(abs(h_current - h_prev))
      
      # Since we normalize h, check relative changes
      if (b_change < 1e-5 && h_change < 1e-5) {
        break
      }
    }
  }

  # Normalize and align HRF shapes
  result <- normalize_and_align_hrf(h_current, b_current, Phi_recon_matrix,
                                   h_ref_shape_canonical, epsilon_scale,
                                   Y_proj, X_list_proj)

  attr(result$h, "iterations") <- iter_final
  attr(result$h, "objective_trace") <- obj_trace[seq_len(iter_final)]
  list(h = result$h, beta = result$beta)
}

