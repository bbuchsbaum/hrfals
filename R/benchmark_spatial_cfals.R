#' Benchmark Spatial CF-ALS
#'
#' Runs the spatial CF-ALS solver for a grid of `lambda_s` values and
#' solver choices on a small synthetic dataset. The function returns a
#' data frame summarising runtime and estimation error for each
#' configuration. It is intended for quick performance exploration as
#' described in HALS-S08 of the proposal.
#'
#' @param lambda_s_values Numeric vector of spatial penalty strengths to test.
#' @param solver_options Character vector of solver options ("direct" or "cg").
#' @param v_seq Vector of voxel counts for the synthetic dataset.
#' @param d Number of basis functions.
#' @param n Number of time points per voxel.
#' @return Data frame with columns `v`, `lambda_s`, `solver`, `runtime`,
#'   and `h_rmse`.
#' @export
benchmark_spatial_cfals <- function(lambda_s_values = c(0, 0.05, 0.1),
                                    solver_options = c("direct", "cg"),
                                    v_seq = c(20, 40),
                                    d = 3,
                                    n = 60) {
  k <- 2
  results <- data.frame()
  for (v in v_seq) {
    set.seed(123 + v)
    h_true <- matrix(rnorm(d * v), d, v)
    beta_true <- matrix(rnorm(k * v), k, v)
    X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
    Y <- matrix(0, n, v)
    for (c in seq_len(k)) {
      Y <- Y + (X_list[[c]] %*% h_true) *
        matrix(rep(beta_true[c, ], each = n), n, v)
    }
    Y <- Y + matrix(rnorm(n * v, sd = 0.1), n, v)
    Phi <- diag(d)
    href <- rep(1, d)
    mask <- array(1, dim = c(v, 1, 1))
    lap <- build_voxel_laplacian(mask)
    for (lambda_s in lambda_s_values) {
      for (solver in solver_options) {
        t <- system.time({
          res <- cf_als_engine(X_list, Y,
                                lambda_b = 0.1,
                                lambda_h = 0.1,
                                lambda_s = lambda_s,
                                laplacian_obj = lap,
                                h_solver = solver,
                                cg_max_iter = 50,
                                cg_tol = 1e-4,
                                precompute_xty_flag = TRUE,
                                Phi_recon_matrix = Phi,
                                h_ref_shape_canonical = href,
                                max_alt = 1)
        })["elapsed"]
        h_err <- sqrt(mean((res$h - h_true)^2))
        results <- rbind(results,
                         data.frame(v = v,
                                    lambda_s = lambda_s,
                                    solver = solver,
                                    runtime = as.numeric(t),
                                    h_rmse = h_err))
      }
    }
  }
  rownames(results) <- NULL
  results
}
