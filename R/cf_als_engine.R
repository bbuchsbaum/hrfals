#' Compute LHS Blocks for the h-update
#'
#' Internal helper used by `cf_als_engine` to construct the list of
#' `d x d` left-hand-side matrices for each voxel's h-update step.
#'
#' @param XtX_list list of `k` cross-product matrices `X_l^T X_l`.
#' @param XtX_full_list optional `k x k` list matrix of full cross-products
#'   `X_l^T X_m` when `fullXtX_flag` is `TRUE`.
#' @param XtX_diag optional `k x d` matrix with column-wise sum-of-squares for
#'   each design block. Used when `fullXtX_flag = FALSE` to form diagonal
#'   Hessian entries efficiently.
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
make_lhs_block_list <- function(XtX_list, XtX_full_list, XtX_diag,
                                b_current, lambda_h, lambda_joint,
                                R_mat_eff, fullXtX_flag, d, v, k) {
  penalty_mat <- if (is.null(R_mat_eff)) diag(d) else R_mat_eff
  base_lhs <- lambda_h * penalty_mat + lambda_joint * diag(d)
  lhs_list <- vector("list", v)

  if (fullXtX_flag) {
    XtX_full_vec <- as.vector(XtX_full_list)
  }

  for (vx in seq_len(v)) {
    lhs_vx <- base_lhs
    b_vx <- b_current[, vx]

    if (fullXtX_flag) {
      weights <- as.vector(outer(b_vx, b_vx))
      weighted_terms <- Map(`*`, XtX_full_vec, weights)
      lhs_vx <- lhs_vx + Reduce(`+`, weighted_terms)
    } else {
      if (!is.null(XtX_diag)) {
        diag(lhs_vx) <- diag(lhs_vx) + drop((b_vx^2) %*% XtX_diag)
      } else {
        weighted_terms <- Map(function(X, w) X * w, XtX_list, b_vx^2)
        lhs_vx <- lhs_vx + Reduce(`+`, weighted_terms)
      }
    }

    lhs_list[[vx]] <- lhs_vx
  }

  lhs_list
}

#' Confound-Adjusted ALS HRF Estimation Engine
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
#' @importFrom stats coef
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

  # Effective beta penalty: when a sparse penalty is used via `beta_penalty$l1`
  # the `lam_beta` ridge component is ignored as described in ticket
  # SP-CFALS-006.
  lambda_b_eff <- if (beta_penalty$l1 > 0) 0 else lambda_b

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
  if (!fullXtX_flag) {
    XtX_diag <- t(vapply(X_list_proj, function(X) colSums(X^2), numeric(d)))
  } else {
    XtX_diag <- NULL
  }


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

    XtY_iter_list <- if (isTRUE(precompute_xty_flag)) {
      XtY_list
    } else {
      lapply(X_list_proj, function(X) crossprod(X, Y_proj))
    }

    DhTy_mat <- matrix(0.0, k, v)
    for (c in seq_len(k)) {
      DhTy_mat[c, ] <- colSums(h_current * XtY_iter_list[[c]])
    }

    G_array <- array(0.0, dim = c(k, k, v))
    if (fullXtX_flag) {
      for (l in seq_len(k)) {
        for (m in seq_len(k)) {
          G_array[l, m, ] <- colSums(h_current *
                                      (XtX_full_list[[l, m]] %*% h_current))
        }
      }
    } else {
      for (l in seq_len(k)) {
        G_array[l, l, ] <- colSums(h_current * (XtX_list[[l]] %*% h_current))
      }
    }

    if (beta_penalty$l1 > 0) {
      Xh_list <- lapply(X_list_proj, function(X) X %*% h_current)
      if (perform_warm_start && iter == 1) {
        warm_beta <- matrix(0.0, k, v)
        for (vx in seq_len(v)) {
          G_vx <- G_array[, , vx] + lambda_joint * diag(k)
          warm_beta[, vx] <- cholSolve(G_vx + lambda_b_eff * diag(k),
                                       DhTy_mat[, vx],
                                       eps = max(epsilon_svd, epsilon_scale))
        }
        b_current <- warm_beta
      }
    }

    for (vx in seq_len(v)) {
      h_vx <- h_current[, vx]

	      if (beta_penalty$l1 > 0) {
	        Xh_mat <- vapply(seq_len(k), function(c) Xh_list[[c]][, vx], numeric(n))
	        y_vx <- Y_proj[, vx]
	        fit <- tryCatch(
	          glmnet::glmnet(
            x = Xh_mat,
            y = y_vx,
            alpha = beta_penalty$alpha,
            lambda = beta_penalty$l1,
            standardize = FALSE,
            intercept = FALSE
          ),
          error = function(e) e
        )
        if (inherits(fit, "error")) {
          warning(sprintf(
            "glmnet failed for voxel %d; falling back to ridge beta update (%s).",
            vx, conditionMessage(fit)
          ), call. = FALSE)
          DhTy_vx <- DhTy_mat[, vx]
          G_vx <- G_array[, , vx] + lambda_joint * diag(k)
          ridge <- max(lambda_b, 0) + max(lambda_joint, 0)
          b_current[, vx] <- cholSolve(G_vx + max(ridge, epsilon_scale) * diag(k),
                                       DhTy_vx,
                                       eps = max(epsilon_svd, epsilon_scale))
        } else {
          beta_hat <- tryCatch(coef(fit, s = beta_penalty$l1),
                               error = function(e) e)
	          if (inherits(beta_hat, "error")) {
	            warning(sprintf(
	              "glmnet coef() failed for voxel %d; falling back to ridge beta update (%s).",
	              vx, conditionMessage(beta_hat)
	            ), call. = FALSE)
	            DhTy_vx <- DhTy_mat[, vx]
	            G_vx <- G_array[, , vx] + lambda_joint * diag(k)
	            ridge <- max(lambda_b, 0) + max(lambda_joint, 0)
	            b_current[, vx] <- cholSolve(G_vx + max(ridge, epsilon_scale) * diag(k),
	                                         DhTy_vx,
	                                         eps = max(epsilon_svd, epsilon_scale))
	          } else {
	            beta_hat_vec <- as.numeric(beta_hat)
	            beta_hat_names <- rownames(beta_hat)
	            if (!is.null(beta_hat_names) && length(beta_hat_names) == length(beta_hat_vec) &&
	                identical(beta_hat_names[[1]], "(Intercept)")) {
	              beta_hat_vec <- beta_hat_vec[-1]
	            } else if (length(beta_hat_vec) == (k + 1)) {
	              beta_hat_vec <- beta_hat_vec[-1]
	            }
	            if (length(beta_hat_vec) != k || anyNA(beta_hat_vec) || any(!is.finite(beta_hat_vec))) {
	              warning(sprintf(
	                "glmnet returned invalid coefficient vector for voxel %d; falling back to ridge beta update.",
	                vx
	              ), call. = FALSE)
	              DhTy_vx <- DhTy_mat[, vx]
	              G_vx <- G_array[, , vx] + lambda_joint * diag(k)
	              ridge <- max(lambda_b, 0) + max(lambda_joint, 0)
	              b_current[, vx] <- cholSolve(G_vx + max(ridge, epsilon_scale) * diag(k),
	                                           DhTy_vx,
	                                           eps = max(epsilon_svd, epsilon_scale))
	            } else {
	              b_current[, vx] <- beta_hat_vec
	            }
	          }
	        }
	      } else {
	        DhTy_vx <- DhTy_mat[, vx]
	        G_vx <- G_array[, , vx] + lambda_joint * diag(k)
        b_current[, vx] <- cholSolve(G_vx + lambda_b_eff * diag(k), DhTy_vx,
                                     eps = max(epsilon_svd, epsilon_scale))
      }
    }

    lhs_block_list <- make_lhs_block_list(
      XtX_list, XtX_full_list, XtX_diag,
      b_current, lambda_h, lambda_joint, R_mat_eff,
      fullXtX_flag, d, v, k
    )

    if (lambda_s > 0 && current_solver == "cg") {
      RHS_mat <- matrix(0.0, d, v)
      for (l in seq_len(k)) {
        RHS_mat <- RHS_mat + XtY_iter_list[[l]] *
          matrix(rep(b_current[l, ], each = d), d, v)
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
      for (l in seq_len(k)) {
        RHS_mat <- RHS_mat + XtY_iter_list[[l]] *
          matrix(rep(b_current[l, ], each = d), d, v)
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
      RHS_mat <- matrix(0.0, d, v)
      for (l in seq_len(k)) {
        RHS_mat <- RHS_mat + XtY_iter_list[[l]] *
          matrix(rep(b_current[l, ], each = d), d, v)
      }
      for (vx in seq_len(v)) {
        lhs <- lhs_block_list[[vx]]
        rhs <- RHS_mat[, vx]
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
    beta_pen <- lambda_b_eff * sum(b_current^2)
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
