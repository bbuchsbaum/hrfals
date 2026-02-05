#' Benchmark fastLSS implementation
#'
#' Provides simple performance benchmarks comparing the fastLSS
#' implementation in \code{fastlss_shared()} to a naive R reference.
#' The function also supports basic scaling tests over different
#' dataset sizes and reports approximate memory usage.
#'
#' @param n_seq Vector of time points to test.
#' @param T_seq Vector of trial counts to test.
#' @param v_seq Vector of voxel counts to test.
#' @param lambda_ridge Optional ridge penalty applied in both
#'   implementations.
#' @return A data frame summarising runtimes and speed ups.
#' @keywords internal
benchmark_fastlss <- function(n_seq = c(40, 80),
                              T_seq = c(10, 20),
                              v_seq = c(5, 10),
                              lambda_ridge = 0) {
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

  simulate_data <- function(n, Tt, v, m = 3) {
    set.seed(1)
    A <- matrix(stats::rnorm(n * m), n, m)
    C <- matrix(stats::rnorm(n * Tt), n, Tt)
    Y <- matrix(stats::rnorm(n * v), n, v)
    p_vec <- stats::rnorm(n)
    list(Y = Y, A = A, C = C, p = p_vec)
  }

  results <- data.frame()
  for (n in n_seq) {
    for (Tt in T_seq) {
      for (v in v_seq) {
        dat <- simulate_data(n, Tt, v)
        mem <- sum(sapply(dat, utils::object.size))
        t_fast <- system.time(
          fastlss_shared(dat$Y, dat$A, dat$C, dat$p,
                         lambda_ridge = lambda_ridge)
        )["elapsed"]
        t_naive <- system.time(
          naive_lss_mode_a(dat$Y, dat$A, dat$C, dat$p,
                            lambda_ridge = lambda_ridge)
        )["elapsed"]
        results <- rbind(results,
                         data.frame(n = n, T = Tt, v = v,
                                    fast = as.numeric(t_fast),
                                    naive = as.numeric(t_naive),
                                    memory = as.numeric(mem)))
      }
    }
  }
  results$speedup <- results$naive / results$fast
  results
}
