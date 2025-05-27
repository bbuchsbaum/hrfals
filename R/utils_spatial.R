#' Build a Voxel Graph Laplacian
#'
#' Constructs an unnormalised graph Laplacian matrix for the voxels
#' contained in a 3D mask volume. The function accepts objects that can
#' be coerced to an array (e.g. a `neuroim2::NeuroVol`). Non-zero
#' elements of the volume are treated as included voxels.
#'
#' @param volume 3D array or `neuroim2::NeuroVol` mask defining the
#'   spatial domain.
#' @param connectivity Neighbourhood definition. Either 6, 18 or 26.
#' @return List with sparse Laplacian matrix `L` (`dgCMatrix`) and the
#'   numeric degree vector `degree`.
#' @export
build_voxel_laplacian <- function(volume, connectivity = 6) {
  build_voxel_laplacian_cpp(volume, as.integer(connectivity))
}

#' Construct Spatial h-update System Matrix
#'
#' Helper used by the spatial CF-ALS implementation. Forms the sparse
#' block matrix \eqn{A_H = BlockDiag(LHS_v) + \lambda_s (L \otimes I_d)}.
#'
#' @param lhs_block_list list of length `v` with `d x d` matrices for each voxel.
#' @param lambda_s spatial regularization strength.
#' @param L_mat sparse Laplacian matrix (`v x v`).
#' @param d number of basis functions.
#' @param v number of voxels.
#' @return `dgCMatrix` representing `A_H`.
#' @keywords internal
construct_A_H_sparse <- function(lhs_block_list, lambda_s, L_mat, d, v) {
  A_block <- Matrix::bdiag(lhs_block_list)
  if (lambda_s > 0) {
    kron_part <- lambda_s * Matrix::kronecker(L_mat, Matrix::Diagonal(d))
    A_block <- A_block + kron_part
  }
  A_block
}

#' Construct Block-Jacobi Preconditioner
#'
#' Forms a block-diagonal preconditioner matrix for the conjugate gradient solver.
#'
#' @param lhs_block_list list of per-voxel `d x d` matrices.
#' @param lambda_s spatial regularization strength.
#' @param degree_vec numeric degree vector from the Laplacian.
#' @param d number of basis functions.
#' @param v number of voxels.
#' @return block-diagonal `dgCMatrix` preconditioner.
#' @keywords internal
construct_preconditioner <- function(lhs_block_list, lambda_s, degree_vec, d, v) {
  pre_list <- vector("list", v)
  for (vx in seq_len(v)) {
    pre_list[[vx]] <- solve(lhs_block_list[[vx]] +
                             lambda_s * degree_vec[vx] * diag(d))
  }
  Matrix::bdiag(pre_list)
}

