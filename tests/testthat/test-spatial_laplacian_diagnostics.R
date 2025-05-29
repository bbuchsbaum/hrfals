context("spatial Laplacian diagnostics")

library(Matrix)

test_that("build_voxel_laplacian edge cases work correctly", {
  # Empty mask should error
  vol_empty <- array(0, dim = c(1, 1, 1))
  expect_error(build_voxel_laplacian(vol_empty), "mask contains no voxels")
  
  # 2D input should error
  vol_2d <- matrix(1, 2, 2)
  expect_error(build_voxel_laplacian(vol_2d), "must be 3D")
})

test_that("cf_als_engine parameter validation works", {
  simple_dat <- simple_small_data()
  
  # Negative lambda_s should error
  expect_error(
    cf_als_engine(simple_dat$X_list, simple_dat$Y,
                  lambda_b = 0.1,
                  lambda_h = 0.2,
                  lambda_s = -0.01,
                  Phi_recon_matrix = simple_dat$Phi,
                  h_ref_shape_canonical = simple_dat$href),
    "lambda_s must be non-negative"
  )
  
  # Missing laplacian_obj should error
  expect_error(
    cf_als_engine(simple_dat$X_list, simple_dat$Y,
                  lambda_b = 0.1,
                  lambda_h = 0.2,
                  lambda_s = 0.05,
                  Phi_recon_matrix = simple_dat$Phi,
                  h_ref_shape_canonical = simple_dat$href),
    "laplacian_obj with elements L and degree must be provided"
  )
  
  # Mismatched degree vector length should error
  bad_lap <- build_voxel_laplacian(array(1, dim = c(2, 1, 1)))
  bad_lap$degree <- bad_lap$degree[1]
  
  expect_error(
    cf_als_engine(simple_dat$X_list, simple_dat$Y,
                  lambda_b = 0.1,
                  lambda_h = 0.2,
                  lambda_s = 0.05,
                  laplacian_obj = bad_lap,
                  Phi_recon_matrix = simple_dat$Phi,
                  h_ref_shape_canonical = simple_dat$href),
    "degree length mismatch"
  )
})

test_that("construct_A_H_sparse and construct_preconditioner work correctly", {
  # Test sparse matrix construction
  L_small <- Matrix(matrix(c(1, -1, -1, 1), 2, 2), sparse = TRUE)
  lhs_list <- list(matrix(2), matrix(2))
  A_res <- hrfals:::construct_A_H_sparse(lhs_list, 0.5, L_small, 1, 2)
  
  expected_A <- Matrix(matrix(c(2.5, -0.5, -0.5, 2.5), 2, 2), sparse = TRUE)
  expect_equal(as.matrix(A_res), as.matrix(expected_A))
  
  # Test preconditioner construction
  P_res <- hrfals:::construct_preconditioner(lhs_list, 0.5, c(1, 1), 1, 2)
  expected_P <- Matrix::bdiag(matrix(1/(2.5), 1, 1), matrix(1/(2.5), 1, 1))
  expect_equal(as.matrix(P_res), as.matrix(expected_P))
})

