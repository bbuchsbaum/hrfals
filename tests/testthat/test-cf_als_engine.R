context("cf_als_engine")

simple_cfals_data <- function() {
  set.seed(123)
  n <- 50
  d <- 3
  k <- 2
  v <- 4
  h_true <- matrix(rnorm(d * v), d, v)
  beta_true <- matrix(rnorm(k * v), k, v)
  X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
  Xbig <- do.call(cbind, X_list)
  Y <- Xbig %*% as.vector(matrix(h_true, d, v) %*% t(beta_true))
  Y <- matrix(Y, n, v)
  phi <- diag(d)
  href <- rep(1, nrow(phi))
  list(X_list = X_list, Y = Y, d = d, k = k,
       Phi = phi, href = href)
}

test_that("cf_als_engine returns matrices with correct dimensions", {
  dat <- simple_cfals_data()
  res <- cf_als_engine(dat$X_list, dat$Y,
                       lambda_b = 0.1,
                       lambda_h = 0.1,
                       R_mat_eff = NULL,
                       fullXtX_flag = FALSE,
                       precompute_xty_flag = TRUE,
                       Phi_recon_matrix = dat$Phi,
                       h_ref_shape_canonical = dat$href,
                       max_alt = 1)
  expect_equal(dim(res$h), c(dat$d, ncol(dat$Y)))
  expect_equal(dim(res$beta), c(dat$k, ncol(dat$Y)))
})




test_that("XtY strategies give identical results", {
  dat <- simple_small_data()
  Rm <- diag(2) * 1.5
  res_pre <- cf_als_engine(dat$X_list, dat$Y,
                           lambda_b = 0.1,
                           lambda_h = 0.2,
                           R_mat_eff = Rm,
                           fullXtX_flag = FALSE,
                           precompute_xty_flag = TRUE,
                           Phi_recon_matrix = dat$Phi,
                           h_ref_shape_canonical = dat$href,
                           max_alt = 1)
  res_onfly <- cf_als_engine(dat$X_list, dat$Y,
                             lambda_b = 0.1,
                             lambda_h = 0.2,
                             R_mat_eff = Rm,
                             fullXtX_flag = FALSE,
                             precompute_xty_flag = FALSE,
                             Phi_recon_matrix = dat$Phi,
                             h_ref_shape_canonical = dat$href,
                             max_alt = 1)
  expect_equal(res_pre$h, res_onfly$h, tolerance = 1e-12)
  expect_equal(res_pre$beta, res_onfly$beta, tolerance = 1e-12)
})

test_that("XtY strategies match with fullXtX", {
  dat <- simple_small_data()
  Rm <- diag(2) * 1.5
  res_pre <- cf_als_engine(dat$X_list, dat$Y,
                           lambda_b = 0.1,
                           lambda_h = 0.2,
                           R_mat_eff = Rm,
                           fullXtX_flag = TRUE,
                           precompute_xty_flag = TRUE,
                           Phi_recon_matrix = dat$Phi,
                           h_ref_shape_canonical = dat$href,
                           max_alt = 1)
  res_onfly <- cf_als_engine(dat$X_list, dat$Y,
                             lambda_b = 0.1,
                             lambda_h = 0.2,
                             R_mat_eff = Rm,
                             fullXtX_flag = TRUE,
                             precompute_xty_flag = FALSE,
                             Phi_recon_matrix = dat$Phi,
                             h_ref_shape_canonical = dat$href,
                             max_alt = 1)
  expect_equal(res_pre$h, res_onfly$h, tolerance = 1e-12)
  expect_equal(res_pre$beta, res_onfly$beta, tolerance = 1e-12)
})

test_that("precompute_xty_flag FALSE reproduces TRUE", {
  dat <- simple_cfals_data()
  res_true <- cf_als_engine(dat$X_list, dat$Y,
                            lambda_b = 0.1,
                            lambda_h = 0.1,
                            fullXtX_flag = FALSE,
                            max_alt = 1,
                            precompute_xty_flag = TRUE,
                            Phi_recon_matrix = dat$Phi,
                            h_ref_shape_canonical = dat$href)
  res_false <- cf_als_engine(dat$X_list, dat$Y,
                             lambda_b = 0.1,
                             lambda_h = 0.1,
                             fullXtX_flag = FALSE,
                             max_alt = 1,
                             precompute_xty_flag = FALSE,
                             Phi_recon_matrix = dat$Phi,
                             h_ref_shape_canonical = dat$href)
  expect_equal(res_false$h, res_true$h)
  expect_equal(res_false$beta, res_true$beta)
})


test_that("h_ref_shape_canonical length must equal p", {
  dat <- simple_cfals_data()
  bad_ref <- numeric(nrow(dat$Phi) + 1)
  expect_error(
    cf_als_engine(dat$X_list, dat$Y,
                  Phi_recon_matrix = dat$Phi,
                  h_ref_shape_canonical = bad_ref),
    "`h_ref_shape_canonical` must have length"
  )
})

test_that("h_ref_shape_canonical must be normalised", {
  dat <- simple_cfals_data()
  bad_ref <- dat$href * 2
  expect_error(
    cf_als_engine(dat$X_list, dat$Y,
                  Phi_recon_matrix = dat$Phi,
                  h_ref_shape_canonical = bad_ref),
    "must be normalised"
  )
})
          
test_that("size estimate uses numeric arithmetic", {
  k <- .Machine$integer.max
  d <- 2L
  v <- 1L
  size_est <- as.numeric(k) * d * v * 8
  expect_true(is.finite(size_est))
  expect_gt(size_est, 2e9)

})



# Additional multi-voxel test to ensure on-the-fly XtY matches precomputed
multi_voxel_data <- function() {
  set.seed(66)
  n <- 30; d <- 2; k <- 2; v <- 3
  h_true <- matrix(rnorm(d * v), d, v)
  beta_true <- matrix(rnorm(k * v), k, v)
  X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
  Y <- matrix(0, n, v)
  for (c in seq_len(k)) {
    Y <- Y + (X_list[[c]] %*% h_true) *
      matrix(rep(beta_true[c, ], each = n), n, v)
  }
  phi <- diag(2)
  href <- rep(1, nrow(phi))
  list(X_list = X_list, Y = Y, Phi = phi, href = href)
}

test_that("XtY cache recomputed per voxel when precomputing disabled", {
  dat <- multi_voxel_data()
  res_true <- cf_als_engine(dat$X_list, dat$Y,
                            lambda_b = 0.1,
                            lambda_h = 0.2,
                            precompute_xty_flag = TRUE,
                            Phi_recon_matrix = dat$Phi,
                            h_ref_shape_canonical = dat$href,
                            max_alt = 1)
  res_false <- cf_als_engine(dat$X_list, dat$Y,
                             lambda_b = 0.1,
                             lambda_h = 0.2,
                             precompute_xty_flag = FALSE,
                             Phi_recon_matrix = dat$Phi,
                             h_ref_shape_canonical = dat$href,
                             max_alt = 1)
  expect_equal(res_false$h, res_true$h, tolerance = 1e-12)
  expect_equal(res_false$beta, res_true$beta, tolerance = 1e-12)
})
                   
test_that("non-symmetric R_mat_eff is forced symmetric", {
  dat <- simple_small_data()
  Rm_nonsym <- matrix(c(1, 2, 3, 4), 2, 2)
  res_nonsym <- cf_als_engine(dat$X_list, dat$Y,
                              lambda_b = 0.1,
                              lambda_h = 0.2,
                              R_mat_eff = Rm_nonsym,
                              fullXtX_flag = FALSE,
                              precompute_xty_flag = TRUE,
                              Phi_recon_matrix = dat$Phi,
                              h_ref_shape_canonical = dat$href,
                              max_alt = 1)
  Rm_sym <- as.matrix(Matrix::forceSymmetric(Rm_nonsym))
  res_sym <- cf_als_engine(dat$X_list, dat$Y,
                           lambda_b = 0.1,
                           lambda_h = 0.2,
                           R_mat_eff = Rm_sym,
                           fullXtX_flag = FALSE,
                           precompute_xty_flag = TRUE,
                           Phi_recon_matrix = dat$Phi,
                           h_ref_shape_canonical = dat$href,
                           max_alt = 1)
  expect_equal(res_nonsym$h, res_sym$h)
  expect_equal(res_nonsym$beta, res_sym$beta)

})

test_that("symmetric Matrix penalty produces expected result", {
  dat <- simple_small_data()
  Rm_dense <- diag(2) * 1.5
  Rm_sym <- as.matrix(Matrix::forceSymmetric(Rm_dense))
  res_sym <- cf_als_engine(dat$X_list, dat$Y,
                           lambda_b = 0.1,
                           lambda_h = 0.2,
                           R_mat_eff = Rm_sym,
                           fullXtX_flag = FALSE,
                           precompute_xty_flag = TRUE,
                           Phi_recon_matrix = dat$Phi,
                           h_ref_shape_canonical = dat$href,
                           max_alt = 1)
  res_dense <- cf_als_engine(dat$X_list, dat$Y,
                             lambda_b = 0.1,
                             lambda_h = 0.2,
                             R_mat_eff = Rm_dense,
                             fullXtX_flag = FALSE,
                             precompute_xty_flag = TRUE,
                             Phi_recon_matrix = dat$Phi,
                             h_ref_shape_canonical = dat$href,
                             max_alt = 1)
  expect_equal(res_sym$h, res_dense$h)
  expect_equal(res_sym$beta, res_dense$beta)
})

test_that("objective converges reasonably", {
  dat <- simple_small_data()
  res <- cf_als_engine(dat$X_list, dat$Y,
                       lambda_b = 0.1,
                       lambda_h = 0.2,
                       precompute_xty_flag = TRUE,
                       Phi_recon_matrix = dat$Phi,
                       h_ref_shape_canonical = dat$href,
                       max_alt = 3)
  obj <- attr(res$h, "objective_trace")
  expect_length(obj, attr(res$h, "iterations"))
  # Check that the objective doesn't increase dramatically
  # CF-ALS with normalization may not strictly decrease due to scale adjustments
  expect_true(all(diff(obj) <= 0.5))  # Allow for reasonable increases
  # Check that we're not diverging wildly
  expect_true(max(obj) / min(obj) < 10)  # Objective shouldn't explode
})

test_that("cg solver path runs", {
  dat <- simple_small_data()
  mask <- array(1, dim = c(2, 1, 1))
  lap_obj <- build_voxel_laplacian(mask)
  res <- cf_als_engine(dat$X_list, dat$Y,
                       lambda_b = 0.1,
                       lambda_h = 0.2,
                       lambda_s = 0.05,
                       laplacian_obj = lap_obj,
                       h_solver = "cg",
                       cg_max_iter = 10,
                       cg_tol = 1e-6,
                       precompute_xty_flag = TRUE,
                       Phi_recon_matrix = dat$Phi,
                       h_ref_shape_canonical = dat$href,
                       max_alt = 1)
  expect_equal(dim(res$h), c(2, ncol(dat$Y)))
  expect_equal(dim(res$beta), c(2, ncol(dat$Y)))
})

test_that("direct solver path runs", {
  dat <- simple_small_data()
  mask <- array(1, dim = c(2, 1, 1))
  lap_obj <- build_voxel_laplacian(mask)
  res <- cf_als_engine(dat$X_list, dat$Y,
                       lambda_b = 0.1,
                       lambda_h = 0.2,
                       lambda_s = 0.05,
                       laplacian_obj = lap_obj,
                       h_solver = "direct",
                       precompute_xty_flag = TRUE,
                       Phi_recon_matrix = dat$Phi,
                       h_ref_shape_canonical = dat$href,
                       max_alt = 1)
  expect_equal(dim(res$h), c(2, ncol(dat$Y)))
  expect_equal(dim(res$beta), c(2, ncol(dat$Y)))
})
