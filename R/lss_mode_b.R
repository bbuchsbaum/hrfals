#' Fast LSS Mode B (voxel-specific trial regressors)
#'
#' Implements the voxel-specific HRF variant of the fast least-squares
#' separate (LSS) algorithm described in `raw-data/FastLSS_proposal.md`.
#' Trial regressors for each voxel are constructed from a list of onset
#' matrices and a matrix of HRF basis coefficients. Computation is
#' performed voxel by voxel using BLAS-optimised operations.
#'
#' @param Y Numeric matrix of BOLD data (n x v).
#' @param A Numeric matrix of nuisance regressors (n x m).
#' @param X_onset_list List of length T containing onset design matrices
#'   (n x d each).
#' @param H_allvoxels Numeric matrix of HRF coefficients (d x v).
#' @param p_vec Numeric vector of length n as described in the proposal.
#' @param lambda_ridge Optional ridge penalty when computing the
#'   pseudoinverse of \code{A}.
#' @return Numeric matrix of trial coefficients (T x v).
#' @export
lss_mode_b <- function(Y, A, X_onset_list, H_allvoxels, p_vec,
                       lambda_ridge = 0) {
  stopifnot(is.matrix(Y), is.matrix(A), is.list(X_onset_list),
            is.matrix(H_allvoxels))
  n <- nrow(Y)
  v <- ncol(Y)
  if (nrow(A) != n)
    stop("Y and A must have the same number of rows")
  if (ncol(H_allvoxels) != v)
    stop("H_allvoxels must have as many columns as Y")
  if (length(p_vec) != n)
    stop("p_vec must have length n")
  Tt <- length(X_onset_list)
  for (i in seq_len(Tt)) {
    if (nrow(X_onset_list[[i]]) != n)
      stop("All onset matrices must have n rows")
  }
  m <- ncol(A)
  AtA <- crossprod(A)
  if (lambda_ridge != 0)
    AtA <- AtA + lambda_ridge * diag(m)
  P <- cholSolve(AtA, t(A))

  trial_names <- names(X_onset_list)
  if (is.null(trial_names))
    trial_names <- paste0("T", seq_len(Tt))
  B <- matrix(0, Tt, v)

  for (vx in seq_len(v)) {
    h_v <- H_allvoxels[, vx]
    C_v <- matrix(0, n, Tt)
    for (t in seq_len(Tt)) {
      C_v[, t] <- X_onset_list[[t]] %*% h_v
    }
    U_v <- P %*% C_v
    V_v <- C_v - A %*% U_v
    pc_row <- drop(crossprod(p_vec, C_v))
    cv_row <- colSums(V_v * V_v)
    alpha_row <- numeric(Tt)
    nz <- cv_row > 0
    alpha_row[nz] <- (1 - pc_row[nz]) / cv_row[nz]
    alpha_row[!nz] <- 0
    S_v <- sweep(V_v, 2, alpha_row, "*")
    S_v <- S_v + p_vec
    B[, vx] <- crossprod(S_v, Y[, vx])
  }
  dimnames(B) <- list(trial_names, colnames(Y))
  B
}

