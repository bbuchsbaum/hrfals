context("lss_mode_a")

naive_lss_mode_a <- function(Y, A, C, p_vec, lambda_ridge = 0) {
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

simple_lss_data <- function() {
  set.seed(123)
  n <- 40; m <- 3; Tt <- 5; v <- 4
  A <- matrix(rnorm(n * m), n, m)
  C <- matrix(rnorm(n * Tt), n, Tt)
  Y <- matrix(rnorm(n * v), n, v)
  p_vec <- rnorm(n)
  list(Y = Y, A = A, C = C, p = p_vec)
}

test_that("lss_mode_a matches naive implementation", {
  dat <- simple_lss_data()
  res_fast <- lss_mode_a(dat$Y, dat$A, dat$C, dat$p, lambda_ridge = 0.1)
  res_naive <- naive_lss_mode_a(dat$Y, dat$A, dat$C, dat$p, lambda_ridge = 0.1)
  expect_equal(res_fast, res_naive, tolerance = 1e-12)
})

test_that("lss_mode_a handles collinear trial", {
  dat <- simple_lss_data()
  dat$C[,1] <- dat$A[,1]  # perfectly collinear with A -> zero residual
  res_fast <- lss_mode_a(dat$Y, dat$A, dat$C, dat$p, lambda_ridge = 0)
  expect_true(all(is.finite(res_fast)))
  res_naive <- naive_lss_mode_a(dat$Y, dat$A, dat$C, dat$p, lambda_ridge = 0)
  expect_equal(res_fast, res_naive, tolerance = 1e-12)
})

test_that("lss_mode_a fallback to QR matches naive", {
  dat <- simple_lss_data()
  res_fast <- lss_mode_a(dat$Y, dat$A, dat$C, dat$p,
                         lambda_ridge = 0.1, woodbury_thresh = 1)
  res_naive <- naive_lss_mode_a(dat$Y, dat$A, dat$C, dat$p,
                                lambda_ridge = 0.1)
  expect_equal(res_fast, res_naive, tolerance = 1e-12)
})
