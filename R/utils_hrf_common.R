#' Validate Common HRF Engine Inputs
#'
#' Internal helper to validate inputs common to both LS+SVD engines.
#'
#' @param X_list_proj list of k design matrices (n x d each)
#' @param Y_proj numeric matrix of projected BOLD data (n x v)
#' @param Phi_recon_matrix Reconstruction matrix mapping coefficients to HRF shape (p x d)
#' @param h_ref_shape_canonical Canonical reference HRF shape of length p for sign alignment
#' @return list with dimensions: n, v, d, k
#' @keywords internal
#' @noRd
validate_hrf_engine_inputs <- function(X_list_proj, Y_proj, Phi_recon_matrix, h_ref_shape_canonical) {
  stopifnot(is.list(X_list_proj), length(X_list_proj) >= 1)
  
  n <- nrow(Y_proj)
  v <- ncol(Y_proj)
  d <- ncol(X_list_proj[[1]])
  k <- length(X_list_proj)
  
  # Validate design matrix dimensions
  for (X in X_list_proj) {
    if (nrow(X) != n) stop("Design matrices must have same rows as Y_proj")
    if (ncol(X) != d) stop("All design matrices must have the same column count")
  }
  
  # Validate reconstruction matrix
  if (!is.matrix(Phi_recon_matrix) || ncol(Phi_recon_matrix) != d)
    stop("`Phi_recon_matrix` must be a p x d matrix")
  
  # Validate canonical HRF shape
  if (length(h_ref_shape_canonical) != nrow(Phi_recon_matrix))
    stop("`h_ref_shape_canonical` must have length nrow(Phi_recon_matrix)")
  if (abs(max(abs(h_ref_shape_canonical)) - 1) > 1e-6)
    stop("`h_ref_shape_canonical` must be normalised to have max abs of 1")
  
  list(n = n, v = v, d = d, k = k)
}

#' Normalize and Align HRF Shapes
#'
#' Internal helper to normalize HRF coefficient matrices and align signs
#' with canonical reference shape. This is the final post-processing step
#' common to both LS+SVD engines.
#'
#' @param H_coeff coefficient matrix (d x v)
#' @param B_coeff beta coefficient matrix (k x v)
#' @param Phi_recon_matrix Reconstruction matrix mapping coefficients to HRF shape (p x d)
#' @param h_ref_shape_canonical Canonical reference HRF shape of length p for sign alignment
#' @param epsilon_scale tolerance for scale in identifiability step
#' @param Y_proj original Y matrix for column names
#' @param X_list_proj original X list for row names
#' @return list with normalized matrices h and beta
#' @keywords internal
#' @noRd
normalize_and_align_hrf <- function(H_coeff, B_coeff, Phi_recon_matrix, 
                                   h_ref_shape_canonical, epsilon_scale,
                                   Y_proj, X_list_proj) {
  v <- ncol(H_coeff)
  
  # Reconstruct HRF shapes
  H_shapes <- Phi_recon_matrix %*% H_coeff
  
  # Calculate scales and alignment
  scl <- apply(abs(H_shapes), 2, max)
  flip <- rep(1.0, v)
  align_scores <- colSums(H_shapes * h_ref_shape_canonical)
  flip[align_scores < 0 & scl > epsilon_scale] <- -1.0
  
  # Apply normalization and sign alignment
  eff_scl <- pmax(scl, epsilon_scale)
  H_final <- sweep(H_coeff, 2, flip / eff_scl, "*")
  B_final <- sweep(B_coeff, 2, flip * eff_scl, "*")
  
  # Zero out coefficients for voxels with negligible scale
  zero_idx <- scl <= epsilon_scale
  if (any(zero_idx)) {
    H_final[, zero_idx] <- 0
    B_final[, zero_idx] <- 0
  }
  
  # Set dimension names
  dimnames(H_final) <- list(NULL, colnames(Y_proj))
  dimnames(B_final) <- list(names(X_list_proj), colnames(Y_proj))
  
  list(h = H_final, beta = B_final)
} 