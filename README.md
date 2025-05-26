# hrfals

[![R-CMD-check](https://github.com/bbuchsbaum/hrfals/workflows/R-CMD-check/badge.svg)](https://github.com/bbuchsbaum/hrfals/actions)

## Overview

The `hrfals` package implements the Confound-Free Alternating Least Squares (CF-ALS) method for estimating hemodynamic response functions (HRFs) from fMRI data. This package provides fast, accurate, and robust alternatives for data-driven HRF estimation that work with any HRF basis from the `fmrireg` package.

## Key Features

- **CF-ALS Algorithm**: Implements the rank-1 decomposition model Y ≈ D(h)β^T for simultaneous estimation of HRF coefficients and activation amplitudes
- **Multiple Estimation Methods**: Supports LS+SVD initialization, LS+SVD+1ALS refinement, and full CF-ALS alternating optimization
- **HRF Basis Compatibility**: Works with any HRF basis from `fmrireg` (B-spline, SPM canonical, tent functions, etc.)
- **Confound Projection**: Built-in QR-based orthogonal projection for nuisance regressor removal
- **Regularization**: Separate L2 penalties for amplitude (λ_β) and HRF coefficient (λ_h) updates
- **Efficient Implementation**: Optimized precomputation and vectorized operations

## Installation

You can install the development version of hrfals from GitHub:

```r
# install.packages("devtools")
devtools::install_github("bbuchsbaum/hrfals")
```

## Quick Start

```r
library(hrfals)
library(fmrireg)

# Create sampling frame and event model
sframe <- fmrireg::sampling_frame(blocklens = 40, TR = 1)
ev_df <- data.frame(onset = c(5, 15, 25), block = 1, cond = "A")
emod <- fmrireg::event_model(onset ~ hrf(cond, basis = fmrireg::HRF_SPMG1),
                            data = ev_df, block = ~ block,
                            sampling_frame = sframe)

# Simulate some BOLD data
Y_matrix <- matrix(rnorm(40 * 5), 40, 5) # 40 timepoints, 5 voxels

# Fit using CF-ALS
cfals_fit <- fmrireg_hrf_cfals(
  fmri_data_obj = Y_matrix,
  event_model = emod,
  hrf_basis = fmrireg::HRF_SPMG1,
  lam_beta = 5,
  lam_h = 0.5,
  max_alt = 1
)

print(cfals_fit)
```

## Main Functions

- `fmrireg_hrf_cfals()`: Main user-facing function for CF-ALS HRF estimation
- `fmrireg_cfals()`: General wrapper supporting multiple estimation methods
- `create_cfals_design()`: Design matrix creation leveraging fmrireg functionality
- `estimate_hrf_cfals()`: Lower-level estimation function for single event terms

## Algorithm Details

The CF-ALS method estimates HRFs and activation amplitudes simultaneously using:

1. **Confound Projection**: Removes nuisance regressors using QR decomposition
2. **SVD Initialization**: Robust starting values via regularized least squares + SVD
3. **Alternating Updates**: 
   - β-update: Solves (G + λ_β I)β = D(h)^T y
   - h-update: Solves (LHS + λ_h I)h = RHS
4. **Identifiability**: Normalizes and aligns HRF estimates for consistency

## Dependencies

- R (>= 3.5.0)
- fmrireg (>= 0.0.0.9000)
- stats

## License

GPL-3

## Citation

If you use this package in your research, please cite:

```
Buchsbaum, B. (2024). hrfals: HRF Estimation using Confound-Free Alternating Least Squares. 
R package version 0.0.0.9000. https://github.com/bbuchsbaum/hrfals
```

## Issues and Support

Please report bugs and feature requests at: https://github.com/bbuchsbaum/hrfals/issues 