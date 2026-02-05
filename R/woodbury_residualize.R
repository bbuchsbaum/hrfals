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
#' set.seed(1)
#' n <- 20; m <- 5; T <- 10
#' A <- matrix(rnorm(n * m), n, m)
#' C <- matrix(rnorm(n * T), n, T)
#' V1 <- woodbury_residualize(C, A)
#' # explicit projection for comparison
#' AtA <- crossprod(A)
#' P <- solve(AtA, crossprod(A, C))
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
  
  # Handle case where A has no columns (no confounds to residualize)
  if (m == 0) {
    return(C)
  }
  
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
#' Residualize trial regressors via QR projection
#'
#' When the number of nuisance regressors is large the Woodbury identity
#' offers little advantage. This helper uses QR decomposition to project
#' trial regressors into the null space of \code{A}.
#'
#' @param C Numeric matrix of trial regressors (n \eqn{\times} T).
#' @param A Numeric matrix of nuisance regressors (n \eqn{\times} m).
#' @param lapack_qr Logical; passed to \code{qr()} for improved stability.
#' @return Matrix of residualised regressors.
#' @keywords internal
qr_residualize <- function(C, A, lapack_qr = FALSE) {
  stopifnot(is.matrix(C), is.matrix(A))
  
  # Handle case where A has no columns (no confounds to residualize)
  if (ncol(A) == 0) {
    return(C)
  }
  
  qrA <- qr(A, LAPACK = lapack_qr)
  
  # qr.resid doesn't support LAPACK QR, so use manual projection
  if (lapack_qr) {
    # Manual projection: C - Q * (Q^T * C)
    Q <- qr.Q(qrA)
    C - Q %*% (t(Q) %*% C)
  } else {
    qr.resid(qrA, C)
  }
}

#' Automatically choose residualisation strategy
#'
#' Uses \code{woodbury_residualize()} when \code{ncol(A)} is below the
#' specified threshold and falls back to \code{qr_residualize()} otherwise.
#'
#' @inheritParams woodbury_residualize
#' @param woodbury_thresh Column threshold for switching to the QR method.
#' @param lapack_qr Logical; passed to \code{qr_residualize} when used.
#' @return Matrix of residualised regressors.
#' @export
auto_residualize <- function(C, A, lambda_ridge = 0,
                             woodbury_thresh = 50,
                             lapack_qr = FALSE) {
  if (ncol(A) > woodbury_thresh) {
    qr_residualize(C, A, lapack_qr = lapack_qr)
  } else {
    woodbury_residualize(C, A, lambda_ridge = lambda_ridge)
  }
}

#' Configure BLAS threading
#'
#' Attempts to set the number of BLAS threads using the optional
#' \pkg{RhpcBLASctl} package. Returns the BLAS library name invisibly.
#'
#' @param num_threads Number of threads to use. Defaults to all available cores.
#' @return Invisibly, a character string describing the BLAS library.
#' @export
configure_blas <- function(num_threads = parallel::detectCores()) {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(num_threads)
  }
  invisible(extSoftVersion()["BLAS"])
}
