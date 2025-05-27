#' Residualize trial regressors via Woodbury identity
#'
#' Computes \eqn{V = C - A(A^+C)} where \eqn{A^+ = (A^TA + \lambda I)^{-1}A^T}.
#' This avoids forming the explicit projection matrix and uses BLAS-optimized
#' matrix multiplications.
#'
#' @param C Numeric matrix of trial regressors (n \eqn{\times} T).
#' @param A Numeric matrix of nuisance regressors (n \eqn{\times} m).
#' @param lambda_ridge Ridge penalty for computing the pseudoinverse of
#'   \eqn{A}. Defaults to 0 (no penalty).
#' @return Matrix of the same dimensions as \code{C} containing residualized
#'   regressors.
#' @examples
#' n <- 20; m <- 5; T <- 10
#' A <- matrix(rnorm(n * m), n, m)
#' C <- matrix(rnorm(n * T), n, T)
#' V1 <- woodbury_residualize(C, A)
#' # explicit projection for comparison
#' AtA <- crossprod(A)
#' P <- cholSolve(AtA, crossprod(A, C))
#' V2 <- C - A %*% P
#' stopifnot(max(abs(V1 - V2)) < 1e-12)
#' @export
woodbury_residualize <- function(C, A, lambda_ridge = 0) {
  stopifnot(is.matrix(C), is.matrix(A))
  n <- nrow(A)
  if (nrow(C) != n) {
    stop("C and A must have the same number of rows")
  }
  m <- ncol(A)
  AtA <- crossprod(A)
  if (lambda_ridge != 0) {
    AtA <- AtA + lambda_ridge * diag(m)
  }
  # solve for A^+ C
  U <- cholSolve(AtA, crossprod(A, C))
  V <- C - A %*% U
  dimnames(V) <- dimnames(C)
  V
}
