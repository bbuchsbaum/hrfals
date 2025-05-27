simulate_cfals_wrapper_data <- function(hrf_basis, noise_sd = 0.05, signal_scale = 1) {
  sf <- fmrireg::sampling_frame(blocklens = 60, TR = 1)
  events <- data.frame(
    onset = c(5, 15, 30, 45),
    condition = factor(c("A", "A", "B", "B")),
    block = 1
  )
  emod <- fmrireg::event_model(onset ~ fmrireg::hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  # Use create_fmri_design to properly create design matrices
  design <- create_fmri_design(emod, hrf_basis)
  X_list <- design$X_list
  
  d <- design$d
  k <- design$k
  v <- 2
  n_timepoints <- length(fmrireg::samples(sf, global = TRUE))
  
  h_true <- matrix(rnorm(d * v), d, v) * signal_scale
  beta_true <- matrix(rnorm(k * v), k, v) * signal_scale
  Y <- matrix(0, n_timepoints, v)
  for (c in seq_along(X_list)) {
    Y <- Y + (X_list[[c]] %*% h_true) *
      matrix(rep(beta_true[c, ], each = nrow(Y)), nrow(Y), v)
  }
  Y <- Y + matrix(rnorm(length(Y), sd = noise_sd), nrow(Y), v)
  attr(Y, "sampling_frame") <- sf
  list(Y = Y, event_model = emod, X_list = X_list,
       h_true = h_true, beta_true = beta_true, sframe = sf)
}
