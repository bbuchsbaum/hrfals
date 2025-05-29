context("cf_als_engine block h-update")

test_that("h-update solves block system", {
  dat <- simple_small_data()
  lambda_b <- 0.1
  lambda_h <- 0.2
  lambda_init <- 1
  Rm <- diag(2)

  res_engine <- cf_als_engine(dat$X_list, dat$Y,
                              lambda_b = lambda_b,
                              lambda_h = lambda_h,
                              lambda_init = lambda_init,
                              R_mat_eff = Rm,
                              fullXtX_flag = FALSE,
                              precompute_xty_flag = TRUE,
                              Phi_recon_matrix = dat$Phi,
                              h_ref_shape_canonical = dat$href,
                              max_alt = 1)

  init <- ls_svd_engine(dat$X_list, dat$Y,
                        lambda_init = lambda_init,
                        Phi_recon_matrix = dat$Phi,
                        h_ref_shape_canonical = dat$href)

  h_current <- init$h
  b_current <- init$beta

  XtX_list <- lapply(dat$X_list, crossprod)
  XtY_list <- lapply(dat$X_list, function(X) crossprod(X, dat$Y))
  k <- length(dat$X_list)
  d <- ncol(dat$X_list[[1]])
  v <- ncol(dat$Y)
  lambda_joint <- 0

  for (vx in seq_len(v)) {
    h_vx <- h_current[, vx]
    DhTy_vx <- vapply(seq_len(k), function(c) {
      crossprod(h_vx, XtY_list[[c]][, vx])
    }, numeric(1))
    G_vx <- matrix(0, k, k)
    for (l in seq_len(k)) {
      G_vx[l, l] <- crossprod(h_vx, XtX_list[[l]] %*% h_vx)
    }
    G_vx <- G_vx + lambda_joint * diag(k)
    b_current[, vx] <- cholSolve(G_vx + lambda_b * diag(k), DhTy_vx)
  }

  lhs_block_list <- hrfals:::make_lhs_block_list(
    XtX_list, NULL, b_current,
    lambda_h, lambda_joint, Rm,
    FALSE, d, v, k
  )

  for (vx in seq_len(v)) {
    rhs <- numeric(d)
    for (l in seq_len(k)) {
      rhs <- rhs + b_current[l, vx] * XtY_list[[l]][, vx]
    }
    h_current[, vx] <- cholSolve(lhs_block_list[[vx]], rhs)
    s <- max(abs(h_current[, vx]), 1e-8)
    h_current[, vx] <- h_current[, vx] / s
    b_current[, vx] <- b_current[, vx] * s
  }

  manual <- normalize_and_align_hrf(
    h_current, b_current,
    dat$Phi, dat$href, 1e-8,
    dat$Y, dat$X_list
  )

  expect_equal(res_engine$h, manual$h, tolerance = 1e-12)
  expect_equal(res_engine$beta, manual$beta, tolerance = 1e-12)
})
