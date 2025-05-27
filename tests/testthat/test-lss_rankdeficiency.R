context("lss rank deficiency")

naive_lss_mode_a_rd <- function(Y, A, C, p_vec, lambda_ridge = 0) {
  n <- nrow(Y); m <- ncol(A); Tt <- ncol(C)
  AtA <- crossprod(A)
  if (lambda_ridge != 0) AtA <- AtA + lambda_ridge * diag(m)
  P <- cholSolve(AtA, t(A))
  B <- matrix(0, Tt, ncol(Y))
  for (t in seq_len(Tt)) {
    c_t <- C[, t]
    u_t <- P %*% c_t
    v_t <- c_t - A %*% u_t
    pc <- crossprod(p_vec, c_t)
    cv <- sum(v_t^2)
    alpha <- if (cv > 0) (1 - pc) / cv else 0
    s_t <- p_vec + as.numeric(alpha) * v_t
    B[t, ] <- crossprod(s_t, Y)
  }
  dimnames(B) <- list(colnames(C), colnames(Y))
  B
}

naive_lss_mode_b_rd <- function(Y, A, X_onset_list, H_allvoxels, p_vec,
                                lambda_ridge = 0) {
  n <- nrow(Y); m <- ncol(A); Tt <- length(X_onset_list); v <- ncol(Y)
  AtA <- crossprod(A)
  if (lambda_ridge != 0) AtA <- AtA + lambda_ridge * diag(m)
  P <- cholSolve(AtA, t(A))
  B <- matrix(0, Tt, v)
  for (vx in seq_len(v)) {
    h_v <- H_allvoxels[, vx]
    C_v <- matrix(0, n, Tt)
    for (t in seq_len(Tt)) {
      C_v[, t] <- X_onset_list[[t]] %*% h_v
      c_t <- C_v[, t]
      u_t <- P %*% c_t
      v_t <- c_t - A %*% u_t
      pc <- crossprod(p_vec, c_t)
      cv <- sum(v_t^2)
      alpha <- if (cv > 0) (1 - pc) / cv else 0
      s_t <- p_vec + as.numeric(alpha) * v_t
      B[t, vx] <- crossprod(s_t, Y[, vx])
    }
  }
  dimnames(B) <- list(NULL, colnames(Y))
  B
}

set.seed(123)

test_that("lss_mode_a handles rank-deficient A", {
  n <- 30; m <- 4; Tt <- 3; v <- 2
  A <- matrix(rnorm(n * m), n, m)
  A[,4] <- A[,1]  # introduce linear dependency
  C <- matrix(rnorm(n * Tt), n, Tt)
  Y <- matrix(rnorm(n * v), n, v)
  p_vec <- rnorm(n)
  lambda <- 0.1
  res_fast <- lss_mode_a(Y, A, C, p_vec, lambda_ridge = lambda)
  res_naive <- naive_lss_mode_a_rd(Y, A, C, p_vec, lambda_ridge = lambda)
  expect_equal(res_fast, res_naive, tolerance = 1e-9)
  expect_true(all(is.finite(res_fast)))
})

set.seed(321)

test_that("lss_mode_b handles rank-deficient A", {
  n <- 25; m <- 4; Tt <- 3; v <- 2; d <- 2
  A <- matrix(rnorm(n * m), n, m)
  A[,4] <- A[,2]  # linear dependency
  X_onset_list <- replicate(Tt, matrix(rnorm(n * d), n, d), simplify = FALSE)
  H_allvoxels <- matrix(rnorm(d * v), d, v)
  Y <- matrix(rnorm(n * v), n, v)
  p_vec <- rnorm(n)
  lambda <- 0.1
  res_fast <- lss_mode_b(Y, A, X_onset_list, H_allvoxels, p_vec,
                         lambda_ridge = lambda)
  res_naive <- naive_lss_mode_b_rd(Y, A, X_onset_list, H_allvoxels, p_vec,
                                   lambda_ridge = lambda)
  expect_equal(res_fast, res_naive, tolerance = 1e-9)
  expect_true(all(is.finite(res_fast)))
})
