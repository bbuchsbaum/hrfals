simulate_cfals_wrapper_data <- function(hrf_basis, noise_sd = 0.05, signal_scale = 1) {
  sf <- sampling_frame(blocklens = 60, TR = 1)
  events <- data.frame(
    onset = c(5, 15, 30, 45),
    condition = factor(c("A", "A", "B", "B")),
    block = 1
  )
  emod <- event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  reg_lists <- lapply(emod$terms, regressors.event_term,
                      hrf = hrf_basis,
                      sampling_frame = sf,
                      summate = FALSE,
                      drop.empty = TRUE)
  regs <- unlist(reg_lists, recursive = FALSE)
  sample_times <- samples(sf, global = TRUE)
  X_list <- lapply(regs, function(r)
    evaluate(r, sample_times, precision = sf$precision))
  d <- nbasis(hrf_basis)
  k <- length(X_list)
  v <- 2
  h_true <- matrix(rnorm(d * v), d, v) * signal_scale
  beta_true <- matrix(rnorm(k * v), k, v) * signal_scale
  Y <- matrix(0, nrow(sample_times), v)
  for (c in seq_along(X_list)) {
    Y <- Y + (X_list[[c]] %*% h_true) *
      matrix(rep(beta_true[c, ], each = nrow(Y)), nrow(Y), v)
  }
  Y <- Y + matrix(rnorm(length(Y), sd = noise_sd), nrow(Y), v)
  attr(Y, "sampling_frame") <- sf
  list(Y = Y, event_model = emod, X_list = X_list,
       h_true = h_true, beta_true = beta_true, sframe = sf)
}
