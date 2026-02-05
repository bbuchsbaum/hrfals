#' Solve linear systems via Cholesky with ridge fallback
#'
#' Attempts to solve \code{M \%*\% x = B} using Cholesky factorisation. If
#' \code{chol()} fails or the resulting factor has a very small diagonal,
#' a ridge of size \code{eps} is added to the diagonal of \code{M}. The
#' ridge factor is scaled adaptively and retried a limited number of times
#' before falling back to a QR-based solve.
#'
#' @param M Symmetric positive definite matrix.
#' @param B Right-hand side matrix or vector.
#' @param eps Base ridge factor added if \code{M} is ill-conditioned.
#' @param max_iter Maximum number of adaptive ridge inflation attempts before
#'   using a QR fallback.
#' @return Solution with the same shape as \code{B}.
#' @keywords internal
#' @noRd
cholSolve <- function(M, B, eps = 1e-8, max_iter = 5L) {
  attempt_chol <- function(mat) {
    tryCatch(chol(mat), error = function(e) NULL)
  }

  L <- attempt_chol(M)
  if (!is.null(L) && min(diag(L)) > eps) {
    return(backsolve(L, forwardsolve(t(L), B)))
  }

  diag_vals <- try(diag(M), silent = TRUE)
  scale_factor <- if (inherits(diag_vals, "try-error") ||
                     any(!is.finite(diag_vals))) 1 else max(1, mean(abs(diag_vals)))
  ridge <- eps * scale_factor
  ident <- diag(nrow(M))

  for (i in seq_len(max_iter)) {
    L <- attempt_chol(M + ridge * ident)
    if (!is.null(L) && min(diag(L)) > sqrt(.Machine$double.eps) * scale_factor) {
      return(backsolve(L, forwardsolve(t(L), B)))
    }
    ridge <- ridge * 10
  }

  qr.solve(M + ridge * ident, B)
}
