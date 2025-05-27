context("spatial Laplacian diagnostics")

library(Matrix)

# build_voxel_laplacian edge cases

vol_empty <- array(0, dim = c(1, 1, 1))

expect_error(build_voxel_laplacian(vol_empty), "mask contains no voxels")

vol_2d <- matrix(1, 2, 2)

expect_error(build_voxel_laplacian(vol_2d), "must be 3D")

# cf_als_engine parameter validation
simple_dat <- simple_small_data()

expect_error(
  cf_als_engine(simple_dat$X_list, simple_dat$Y,
                lambda_b = 0.1,
                lambda_h = 0.2,
                lambda_s = -0.01,
                Phi_recon_matrix = simple_dat$Phi,
                h_ref_shape_canonical = simple_dat$href),
  "lambda_s must be non-negative"
)

expect_error(
  cf_als_engine(simple_dat$X_list, simple_dat$Y,
                lambda_b = 0.1,
                lambda_h = 0.2,
                lambda_s = 0.05,
                Phi_recon_matrix = simple_dat$Phi,
                h_ref_shape_canonical = simple_dat$href),
  "laplacian_obj with elements L and degree must be provided"
)

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
  "laplacian_obj\$degree length mismatch"
)

# construct_A_H_sparse and construct_preconditioner behaviour
L_small <- Matrix(matrix(c(1, -1, -1, 1), 2, 2), sparse = TRUE)
lhs_list <- list(matrix(2), matrix(2))
A_res <- hrfals:::construct_A_H_sparse(lhs_list, 0.5, L_small, 1, 2)

expected_A <- Matrix(matrix(c(2.5, -0.5, -0.5, 2.5), 2, 2), sparse = TRUE)
expect_equal(as.matrix(A_res), as.matrix(expected_A))

P_res <- hrfals:::construct_preconditioner(lhs_list, 0.5, c(1, 1), 1, 2)
expected_P <- Matrix::bdiag(matrix(1/(2.5), 1, 1), matrix(1/(2.5), 1, 1))
expect_equal(as.matrix(P_res), as.matrix(expected_P))

