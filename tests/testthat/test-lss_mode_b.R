context("lss_mode_b")

naive_lss_mode_b <- function(Y, A, X_onset_list, H_allvoxels, p_vec,
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

simple_lss_b_data <- function() {
  set.seed(123)
  n <- 30; m <- 3; Tt <- 4; v <- 3; d <- 2
  A <- matrix(rnorm(n * m), n, m)
  X_onset_list <- replicate(Tt, matrix(rnorm(n * d), n, d), simplify = FALSE)
  H_allvoxels <- matrix(rnorm(d * v), d, v)
  Y <- matrix(rnorm(n * v), n, v)
  p_vec <- rnorm(n)
  list(Y=Y,A=A,X=X_onset_list,H=H_allvoxels,p=p_vec)
}

test_that("lss_mode_b matches naive implementation", {
  dat <- simple_lss_b_data()
  res_fast <- lss_mode_b(dat$Y, dat$A, dat$X, dat$H, dat$p, lambda_ridge = 0.1)
  res_naive <- naive_lss_mode_b(dat$Y, dat$A, dat$X, dat$H, dat$p, lambda_ridge = 0.1)
  expect_equal(unname(res_fast), unname(res_naive), tolerance = 1e-12)
})

test_that("lss_mode_b handles collinear trial", {
  dat <- simple_lss_b_data()
  # make first trial's first basis function perfectly collinear with A's first regressor
  dat$X[[1]][,1] <- dat$A[,1]
  # H should remain d x v, where d is the original number of basis functions (2)
  # No change to dat$H here to maintain dimensional consistency with dat$X[[t]] for t>1
  res_fast <- lss_mode_b(dat$Y, dat$A, dat$X, dat$H, dat$p, lambda_ridge = 0)
  res_naive <- naive_lss_mode_b(dat$Y, dat$A, dat$X, dat$H, dat$p, lambda_ridge = 0)
  expect_equal(unname(res_fast), unname(res_naive), tolerance = 1e-12)
  expect_true(all(is.finite(res_fast)))
})

test_that("lss_mode_b fallback to QR matches naive", {
  dat <- simple_lss_b_data()
  res_fast <- lss_mode_b(dat$Y, dat$A, dat$X, dat$H, dat$p,
                         lambda_ridge = 0.1, woodbury_thresh = 1)
  res_naive <- naive_lss_mode_b(dat$Y, dat$A, dat$X, dat$H, dat$p,
                                lambda_ridge = 0.1)
  expect_equal(unname(res_fast), unname(res_naive), tolerance = 1e-2)
})

# Prewhitening support
test_that("lss_mode_b handles whitening matrix", {
  dat <- simple_lss_b_data()
  n <- nrow(dat$Y)
  set.seed(2)
  W <- chol(crossprod(matrix(rnorm(n*n), n, n)))
  res_fast <- lss_mode_b(dat$Y, dat$A, dat$X, dat$H, dat$p,
                         lambda_ridge = 0.1, W = W)
  datW <- list(Y = W %*% dat$Y,
               A = W %*% dat$A,
               X = lapply(dat$X, function(X) W %*% X),
               H = dat$H,
               p = drop(W %*% dat$p))
  res_naive <- naive_lss_mode_b(datW$Y, datW$A, datW$X, datW$H, datW$p,
                                lambda_ridge = 0.1)
  expect_equal(unname(res_fast), unname(res_naive), tolerance = 1e-12)
})
