#' Solve linear systems via Cholesky with ridge fallback
#'
#' Attempts to solve \code{M %*% x = B} using Cholesky factorisation. If
#' \code{chol()} fails or the resulting factor has a very small diagonal,
#' a ridge of size \code{eps} is added to the diagonal of \code{M} before
#' solving.
#'
#' @param M Symmetric positive definite matrix.
#' @param B Right-hand side matrix or vector.
#' @param eps Ridge factor added if \code{M} is ill-conditioned.
#' @return Solution with the same shape as \code{B}.
#' @keywords internal
#' @noRd
cholSolve <- function(M, B, eps = 1e-8) {
  ok <- TRUE
  L <- tryCatch(chol(M), error = function(e) { ok <<- FALSE ; NULL })
  if (!ok || (is.null(L) || min(diag(L)) < eps)) {
    L <- chol(M + eps * diag(nrow(M)))
  }
  backsolve(L, forwardsolve(t(L), B))
}
