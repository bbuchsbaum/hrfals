context("lss_kernel_cpp")

# Naive implementation following the FastLSS proposal algorithm
naive_lss_kernel <- function(C, A, Y, p_vec, lambda_ridge = 0) {
  # Step 1: Compute P = (A^T A + lambda I)^(-1) A^T
  AtA <- crossprod(A)
  if (lambda_ridge != 0) AtA <- AtA + lambda_ridge * diag(ncol(A))
  P <- hrfals:::cholSolve(AtA, t(A))
  
  # Step 2: Compute U = PC
  U <- P %*% C
  
  # Step 3: Compute V = C - AU (residualized C)
  V <- C - A %*% U
  
  # Step 4: Compute pc_row = p_vec^T * C
  pc_row <- crossprod(p_vec, C)
  
  # Step 5: Compute cv_row = colSums(V * V)
  cv_row <- colSums(V * V)
  
  # Step 6: Compute alpha_row = (1 - pc_row) / cv_row
  alpha_row <- numeric(ncol(C))
  nz <- cv_row > 0
  alpha_row[nz] <- (1 - pc_row[nz]) / cv_row[nz]
  alpha_row[!nz] <- 0
  
  # Step 7: Construct S matrix
  S <- matrix(0, nrow(C), ncol(C))
  for (t in seq_len(ncol(C))) {
    S[, t] <- p_vec + alpha_row[t] * V[, t]
  }
  
  # Step 8: Compute B = S^T * Y
  B <- crossprod(S, Y)
  
  return(B)
}

simple_kernel_data <- function() {
  set.seed(1)
  n <- 20; m <- 4; Tt <- 5; v <- 3
  C <- matrix(rnorm(n*Tt), n, Tt)
  A <- matrix(rnorm(n*m), n, m)
  Y <- matrix(rnorm(n*v), n, v)
  p_vec <- rnorm(n)
  list(C=C, A=A, Y=Y, p_vec=p_vec)
}

test_that("lss_kernel_cpp matches naive implementation", {
  dat <- simple_kernel_data()
  res_cpp <- lss_kernel_cpp(dat$C, dat$A, dat$Y, dat$p_vec, lambda_ridge = 0.1)
  res_naive <- naive_lss_kernel(dat$C, dat$A, dat$Y, dat$p_vec, lambda_ridge = 0.1)
  expect_equal(res_cpp, res_naive, tolerance = 1e-12)
})
