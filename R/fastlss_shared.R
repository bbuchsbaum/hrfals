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
#' @param woodbury_thresh Threshold for switching from Woodbury to
#'   QR-based residualisation. See \code{auto_residualize}.
#' @param chunk_size Optional chunk size (number of trials) used to
#'   process \code{C} in blocks. Set automatically when
#'   \code{mem_limit} is supplied.
#' @param progress Logical; display a progress bar when processing in
#'   chunks.
#' @param mem_limit Optional memory limit in megabytes for automatic
#'   chunking.
#' @param W Optional whitening matrix to apply to `Y`, `A` and `C`
#'   before running the kernel.
#' @param use_cpp Logical; use the enhanced C++ implementation when TRUE.
#' @return A numeric matrix of trial coefficients (T x v).
#' @export
lss_mode_a <- function(Y, A, C, p_vec, lambda_ridge = 0,
                       woodbury_thresh = 50,
                       chunk_size = NULL,
                       progress = FALSE,
                       mem_limit = NULL,
                       W = NULL,
                       use_cpp = FALSE) {
  stopifnot(is.matrix(Y), is.matrix(A), is.matrix(C))
  n <- nrow(Y)
  if (nrow(A) != n || nrow(C) != n)
    stop("Y, A and C must have the same number of rows")
  if (length(p_vec) != n)
    stop("p_vec must have length n")

  if (!is.null(W)) {
    if (!is.matrix(W) || nrow(W) != n || ncol(W) != n)
      stop("'W' must be an n x n whitening matrix")
    Y <- W %*% Y
    A <- W %*% A
    C <- W %*% C
    p_vec <- drop(W %*% p_vec)
  }

  m <- ncol(A)
  Tt <- ncol(C)

  # Use enhanced C++ implementation if requested
  if (use_cpp) {
    if (!is.null(chunk_size) || !is.null(mem_limit) || progress) {
      warning("C++ implementation ignores chunk_size, mem_limit, and progress arguments")
    }
    B <- lss_kernel_cpp(C, A, Y, p_vec, lambda_ridge = lambda_ridge, shared_C = TRUE)
    dimnames(B) <- list(colnames(C), colnames(Y))
    return(B)
  }

  if (!is.null(mem_limit) && is.null(chunk_size)) {
    bytes_limit <- mem_limit * 1024^2
    max_chunk <- floor(bytes_limit / (8 * n))
    if (max_chunk > 0 && max_chunk < Tt)
      chunk_size <- max_chunk
  }

  if (is.null(chunk_size) || chunk_size >= Tt) {
    V <- auto_residualize(C, A, lambda_ridge,
                          woodbury_thresh = woodbury_thresh)
    pc_row <- drop(crossprod(p_vec, C))     # length T
    cv_row <- colSums(V * V)
    alpha_row <- numeric(Tt)
    nz <- cv_row > 0
    alpha_row[nz] <- (1 - pc_row[nz]) / cv_row[nz]
    alpha_row[!nz] <- 0
    S <- sweep(V, 2, alpha_row, "*")
    S <- S + p_vec
    B <- crossprod(S, Y)    # T x v
  } else {
    idx_list <- split(seq_len(Tt), ceiling(seq_len(Tt) / chunk_size))
    B <- matrix(0, Tt, ncol(Y))
    if (progress)
      pb <- txtProgressBar(min = 0, max = length(idx_list), style = 3)
    i <- 0
    for (idx in idx_list) {
      i <- i + 1
      C_chunk <- C[, idx, drop = FALSE]
      V_chunk <- auto_residualize(C_chunk, A, lambda_ridge,
                                  woodbury_thresh = woodbury_thresh)
      pc_row <- drop(crossprod(p_vec, C_chunk))
      cv_row <- colSums(V_chunk * V_chunk)
      alpha_row <- numeric(length(idx))
      nz <- cv_row > 0
      alpha_row[nz] <- (1 - pc_row[nz]) / cv_row[nz]
      alpha_row[!nz] <- 0
      S_chunk <- sweep(V_chunk, 2, alpha_row, "*")
      S_chunk <- S_chunk + p_vec
      B[idx, ] <- crossprod(S_chunk, Y)
      if (progress) setTxtProgressBar(pb, i)
    }
    if (progress) close(pb)
  }

  dimnames(B) <- list(colnames(C), colnames(Y))
  B
}
