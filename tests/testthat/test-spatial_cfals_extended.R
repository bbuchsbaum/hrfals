context("spatial cfals extended tests")

simple_spatial_data <- function() {
  dat <- simple_small_data()
  mask <- array(1, dim = c(2, 1, 1))
  lap_obj <- build_voxel_laplacian(mask)
  list(dat = dat, lap = lap_obj)
}

# lambda_s = 0 should match non-spatial results

test_that("lambda_s=0 matches baseline", {
  sdat <- simple_spatial_data()
  dat <- sdat$dat
  lap <- sdat$lap
  base <- cf_als_engine(dat$X_list, dat$Y,
                        lambda_b = 0.1,
                        lambda_h = 0.2,
                        precompute_xty_flag = TRUE,
                        Phi_recon_matrix = dat$Phi,
                        h_ref_shape_canonical = dat$href,
                        max_alt = 1)
  spatial_zero <- cf_als_engine(dat$X_list, dat$Y,
                                lambda_b = 0.1,
                                lambda_h = 0.2,
                                lambda_s = 0,
                                laplacian_obj = lap,
                                h_solver = "direct",
                                precompute_xty_flag = TRUE,
                                Phi_recon_matrix = dat$Phi,
                                h_ref_shape_canonical = dat$href,
                                max_alt = 1)
  expect_equal(spatial_zero$h, base$h)
  expect_equal(spatial_zero$beta, base$beta)
})

# direct vs CG solver should give similar results

test_that("direct and CG solvers agree", {
  sdat <- simple_spatial_data()
  dat <- sdat$dat
  lap <- sdat$lap
  direct <- cf_als_engine(dat$X_list, dat$Y,
                          lambda_b = 0.1,
                          lambda_h = 0.2,
                          lambda_s = 0.05,
                          laplacian_obj = lap,
                          h_solver = "direct",
                          precompute_xty_flag = TRUE,
                          Phi_recon_matrix = dat$Phi,
                          h_ref_shape_canonical = dat$href,
                          max_alt = 1)
  cg <- cf_als_engine(dat$X_list, dat$Y,
                      lambda_b = 0.1,
                      lambda_h = 0.2,
                      lambda_s = 0.05,
                      laplacian_obj = lap,
                      h_solver = "cg",
                      cg_max_iter = 50,
                      cg_tol = 1e-8,
                      precompute_xty_flag = TRUE,
                      Phi_recon_matrix = dat$Phi,
                      h_ref_shape_canonical = dat$href,
                      max_alt = 1)
  expect_equal(cg$h, direct$h, tolerance = 1e-5)
  expect_equal(cg$beta, direct$beta, tolerance = 1e-5)
})

# sign alignment across voxels even if one voxel is flipped

test_that("sign alignment works with flipped voxel", {
  sdat <- simple_spatial_data()
  dat <- sdat$dat
  lap <- sdat$lap
  # flip sign of second voxel's data
  dat$Y[, 2] <- -dat$Y[, 2]
  res <- cf_als_engine(dat$X_list, dat$Y,
                       lambda_b = 0.1,
                       lambda_h = 0.2,
                       lambda_s = 0.05,
                       laplacian_obj = lap,
                       h_solver = "direct",
                       precompute_xty_flag = TRUE,
                       Phi_recon_matrix = dat$Phi,
                       h_ref_shape_canonical = dat$href,
                       max_alt = 1)
  recon <- dat$Phi %*% res$h
  expect_true(all(colSums(recon * dat$href) > 0))
})

