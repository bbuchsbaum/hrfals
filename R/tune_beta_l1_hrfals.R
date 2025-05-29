#' Tune L1 Penalty for Sparse Beta Estimation
#'
#' Helper performing a simple grid search over the L1 penalty used for
#' sparse beta estimation in `hrfals`. A subset of voxels is split into
#' train and test sets. For each value in `l1_grid` the CF-ALS solver is
#' run on the training voxels and prediction error on the test voxels is
#' computed. The best L1 value is returned along with the full results
#' table.
#'
#' @param fmri_data_obj_subset Numeric matrix or `fmri_dataset` containing
#'   a subset of voxels used for tuning.
#' @param event_model An `event_model` describing the stimuli.
#' @param hrf_basis HRF basis object.
#' @param l1_grid Numeric vector of L1 penalty strengths to evaluate.
#' @param alpha_value Elastic net mixing parameter passed to `hrfals`.
#' @param n_outer_iterations_cfals Number of CF-ALS iterations to run
#'   during tuning.
#' @param other_hrfals_args Named list of additional arguments forwarded
#'   to [hrfals()].
#' @param cv_voxel_subset_train_prop Proportion of voxels used for the
#'   training split. Remaining voxels are used for testing.
#' @param seed Optional random seed for reproducible splits.
#' @return Data frame with columns `l1`, `train_mse`, `test_mse` and a
#'   `best` flag. The attribute `best_l1` holds the optimal penalty.
#' @export
#' @examples
#' \dontrun{
#' res <- tune_beta_l1_hrfals(Y_small, emod, HRF_SPMG3)
#' res$best_l1
#' }

tune_beta_l1_hrfals <- function(fmri_data_obj_subset,
                                event_model,
                                hrf_basis,
                                l1_grid = 10^seq(-3, 0, length.out = 7),
                                alpha_value = 1,
                                n_outer_iterations_cfals = 1,
                                other_hrfals_args = list(lam_beta = 0.01,
                                                         lam_h = 0.01),
                                cv_voxel_subset_train_prop = 0.7,
                                seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  Y <- if (inherits(fmri_data_obj_subset, "fmri_dataset")) {
    fmrireg::get_data_matrix(fmri_data_obj_subset)
  } else if (is.matrix(fmri_data_obj_subset)) {
    fmri_data_obj_subset
  } else {
    stop("'fmri_data_obj_subset' must be an 'fmri_dataset' or matrix")
  }

  v <- ncol(Y)
  train_idx <- sample(seq_len(v), size = floor(cv_voxel_subset_train_prop * v))
  test_idx <- setdiff(seq_len(v), train_idx)

  results <- data.frame(l1 = l1_grid,
                        train_mse = NA_real_,
                        test_mse = NA_real_)

  for (i in seq_along(l1_grid)) {
    l1_val <- l1_grid[i]

    args_common <- c(list(fmri_data_obj = Y[, train_idx, drop = FALSE],
                          event_model = event_model,
                          hrf_basis = hrf_basis,
                          beta_penalty = list(l1 = l1_val,
                                               alpha = alpha_value,
                                               warm_start = TRUE),
                          max_alt = n_outer_iterations_cfals),
                     other_hrfals_args)
    fit_train <- do.call(hrfals, args_common)
    train_mse <- mean(fit_train$residuals^2)

    args_common$fmri_data_obj <- Y[, test_idx, drop = FALSE]
    fit_test <- do.call(hrfals, args_common)
    test_mse <- mean(fit_test$residuals^2)

    results$train_mse[i] <- train_mse
    results$test_mse[i] <- test_mse
  }

  best_idx <- which.min(results$test_mse)
  results$best <- FALSE
  results$best[best_idx] <- TRUE
  attr(results, "best_l1") <- results$l1[best_idx]
  results
}
