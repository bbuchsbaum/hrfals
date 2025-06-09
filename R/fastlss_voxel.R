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
#' @param woodbury_thresh Threshold for switching from Woodbury to
#'   QR-based residualisation. See \code{auto_residualize}.
#' @param chunk_size Optional chunk size (number of trials) used to
#'   process per voxel when memory limits are a concern. Set
#'   automatically when \code{mem_limit} is supplied.
#' @param progress Logical; display a progress bar over voxels when
#'   \code{TRUE}.
#' @param mem_limit Optional memory limit in megabytes for automatic
#'   chunking.
#' @param W Optional whitening matrix applied to `Y`, `A` and each
#'   onset matrix before running the kernel.
#' @return Numeric matrix of trial coefficients (T x v).
#' @export
lss_mode_b <- function(Y, A, X_onset_list, H_allvoxels, p_vec,
                       lambda_ridge = 0,
                       woodbury_thresh = 50,
                       chunk_size = NULL,
                       progress = FALSE,
                       mem_limit = NULL,
                       W = NULL) {
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

  if (!is.null(W)) {
    if (!is.matrix(W) || nrow(W) != n || ncol(W) != n)
      stop("'W' must be an n x n whitening matrix")
    Y <- W %*% Y
    A <- W %*% A
    X_onset_list <- lapply(X_onset_list, function(X) W %*% X)
    p_vec <- drop(W %*% p_vec)
  }

  trial_names <- names(X_onset_list)
  if (is.null(trial_names))
    trial_names <- paste0("T", seq_len(Tt))
  B <- matrix(0, Tt, v)

  ## Precompute stacked onset matrices for efficient voxel loop
  X_onset_stack <- do.call(rbind, X_onset_list)

  if (!is.null(mem_limit) && is.null(chunk_size)) {
    bytes_limit <- mem_limit * 1024^2
    max_chunk <- floor(bytes_limit / (8 * n))
    if (max_chunk > 0 && max_chunk < Tt)
      chunk_size <- max_chunk
  }

  if (progress)
    pb <- txtProgressBar(min = 0, max = v, style = 3)

  for (vx in seq_len(v)) {
    h_v <- H_allvoxels[, vx]
    C_v <- matrix(X_onset_stack %*% h_v, nrow = n, ncol = Tt,
                  byrow = FALSE)

    if (is.null(chunk_size) || chunk_size >= Tt) {
      V_v <- auto_residualize(C_v, A, lambda_ridge,
                              woodbury_thresh = woodbury_thresh)
      pc_row <- drop(crossprod(p_vec, C_v))
      cv_row <- colSums(V_v * V_v)
      alpha_row <- numeric(Tt)
      nz <- cv_row > 0
      alpha_row[nz] <- (1 - pc_row[nz]) / cv_row[nz]
      alpha_row[!nz] <- 0
      S_v <- sweep(V_v, 2, alpha_row, "*")
      S_v <- S_v + p_vec
      B[, vx] <- crossprod(S_v, Y[, vx])
    } else {
      idx_list <- split(seq_len(Tt), ceiling(seq_len(Tt) / chunk_size))
      res_col <- numeric(Tt)
      for (idx in idx_list) {
        C_chunk <- C_v[, idx, drop = FALSE]
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
        res_col[idx] <- crossprod(S_chunk, Y[, vx])
      }
      B[, vx] <- res_col
    }
    if (progress) setTxtProgressBar(pb, vx)
  }
  if (progress) close(pb)
  dimnames(B) <- list(trial_names, colnames(Y))
  B
}

#' @rdname lss_mode_b
#' @export
fastlss_voxel <- lss_mode_b

