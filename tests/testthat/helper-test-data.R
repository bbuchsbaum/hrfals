# Helper functions for test data generation

simple_small_data <- function() {
  set.seed(42)
  n <- 20; d <- 2; k <- 2; v <- 2
  h_true <- matrix(rnorm(d * v), d, v)
  beta_true <- matrix(rnorm(k * v), k, v)
  X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
  Y <- matrix(0, n, v)
  for (c in seq_len(k)) {
    Y <- Y + (X_list[[c]] %*% h_true) *
      matrix(rep(beta_true[c, ], each = n), n, v)
  }
  phi <- diag(d)
  href <- rep(1, nrow(phi))
  list(X_list = X_list, Y = Y, Phi = phi, href = href)
} 