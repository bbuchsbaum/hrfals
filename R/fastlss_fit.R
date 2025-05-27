#' fastlss_fit constructor
#'
#' Simple container for beta-series estimates produced by `hrfals_lss()`.
#' Stores the beta matrix together with metadata about the CF-ALS fit and
#' event information used to build the trial regressors.
#'
#' @param betas Numeric matrix of trial coefficients (T x V).
#' @param mode Character string specifying the LSS kernel variant.
#' @param cfals_fit The `hrfals_fit` object used to derive the HRFs.
#' @param events The event table or `event_model` used when building the
#'   trial regressors.
#' @param hrf_basis HRF basis object employed for the regressors.
#' @param call Matched call to `hrfals_lss()`.
#' @param whitening_matrix Optional whitening matrix applied prior to
#'   estimation.
#' @return An object of class `fastlss_fit`.
#' @export
fastlss_fit <- function(betas, mode, cfals_fit, events, hrf_basis,
                        call = NULL, whitening_matrix = NULL) {
  out <- list(
    betas = betas,
    mode = mode,
    cfals_fit = cfals_fit,
    events = events,
    hrf_basis = hrf_basis,
    call = call,
    whitening_matrix = whitening_matrix
  )
  class(out) <- c("fastlss_fit", "list")
  out
}

#' @export
print.fastlss_fit <- function(x, ...) {
  cat("\nfastLSS beta-series\n")
  cat("===================\n")
  cat(sprintf("Trials: %d\n", nrow(x$betas)))
  cat(sprintf("Voxels: %d\n", ncol(x$betas)))
  cat(sprintf("Mode: %s\n", x$mode))
  invisible(x)
}

#' @export
summary.fastlss_fit <- function(object, ...) {
  structure(
    list(trials = nrow(object$betas),
         voxels = ncol(object$betas),
         mode = object$mode),
    class = "summary.fastlss_fit"
  )
}

#' @export
print.summary.fastlss_fit <- function(x, ...) {
  cat("\nSummary of fastLSS fit\n")
  cat("----------------------\n")
  cat(sprintf("Trials: %d\n", x$trials))
  cat(sprintf("Voxels: %d\n", x$voxels))
  cat(sprintf("Mode: %s\n", x$mode))
  invisible(x)
}
