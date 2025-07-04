% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hrfals_methods.R
\name{predict.hrfals_fit}
\alias{predict.hrfals_fit}
\title{Predict method for hrfals_fit objects}
\usage{
\method{predict}{hrfals_fit}(object, newdata = NULL, ...)
}
\arguments{
\item{object}{An `hrfals_fit` object.}

\item{newdata}{Optional new data for prediction. Currently not implemented.}

\item{...}{Additional arguments (currently unused).}
}
\value{
A matrix of predicted values with dimensions (n_timepoints x n_voxels).
  For the original training data, these are the fitted values from the model.
}
\description{
Generates predictions from a fitted hrfals model. When no new data is provided,
returns the fitted values for the original training data.
}
\details{
The predictions are computed as:
\deqn{Y_{pred} = \sum_{c=1}^{k} (X_c \times h) \odot \beta_c}
where \eqn{X_c} are the design matrices, \eqn{h} are the HRF coefficients,
\eqn{\beta_c} are the condition amplitudes, and \eqn{\odot} denotes 
element-wise multiplication across voxels.
}
\examples{
\dontrun{
# Assuming you have a fitted hrfals object
fitted_values <- predict(fit)

# Check fitted values against original data
residuals_manual <- Y - fitted_values
all.equal(residuals_manual, residuals(fit))
}
}
