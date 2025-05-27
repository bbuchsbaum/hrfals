#' Fast LSS Mode A (shared trial regressors)
#'
#' Implements the shared-regressor variant of the fast least-squares
#' separate (LSS) algorithm described in `raw-data/FastLSS_proposal.md`.
#' All heavy computations are carried out using BLAS-optimised matrix
#' operations.
#'
#' @param Y Numeric matrix of BOLD data (n x v).
#' @param A Numeric matrix of nuisance regressors (n x m).
#' @param C Numeric matrix of trial regressors shared across voxels
#'   (n x T).
#' @param p_vec Numeric vector of length n as described in the proposal.
#' @param lambda_ridge Optional ridge penalty when computing the
#'   pseudoinverse of \code{A}.
#' @return A numeric matrix of trial coefficients (T x v).
#' @export
lss_mode_a <- function(Y, A, C, p_vec, lambda_ridge = 0) {
  stopifnot(is.matrix(Y), is.matrix(A), is.matrix(C))
  n <- nrow(Y)
  if (nrow(A) != n || nrow(C) != n)
    stop("Y, A and C must have the same number of rows")
  if (length(p_vec) != n)
    stop("p_vec must have length n")

  m <- ncol(A)
  Tt <- ncol(C)

  AtA <- crossprod(A)
  if (lambda_ridge != 0)
    AtA <- AtA + lambda_ridge * diag(m)
  P <- cholSolve(AtA, t(A))

  U <- P %*% C            # m x T
  V <- C - A %*% U        # n x T

  pc_row <- drop(crossprod(p_vec, C))     # length T
  cv_row <- colSums(V * V)
  alpha_row <- numeric(Tt)
  nz <- cv_row > 0
  alpha_row[nz] <- (1 - pc_row[nz]) / cv_row[nz]
  alpha_row[!nz] <- 0

  S <- sweep(V, 2, alpha_row, "*")
  S <- S + p_vec

  B <- crossprod(S, Y)    # T x v
  dimnames(B) <- list(colnames(C), colnames(Y))
  B
}
