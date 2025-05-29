# Proposal: Sparse Beta Estimation for Shared HRF with Many Continuous Predictors in the `hrfals` Package

## 1. Motivation and Goal

The `hrfals` package currently provides robust and efficient methods (LS+SVD, LS+SVD+1ALS, CF-ALS) for estimating a shared Hemodynamic Response Function (HRF) and associated condition amplitudes from fMRI data. This proposal outlines an extension to these methods, specifically targeting scenarios where the BOLD signal is modeled by a potentially large number ($K \gg 1$) of continuous predictor time series ($x_k(t)$), all assumed to be convolved with an HRF shape $h_v(\tau)$ that is shared by all $K$ predictors *within a given voxel $v$*. The resulting amplitudes $\beta_{k,v}$ then scale these convolved predictors.

The key innovation is the introduction of sparsity-inducing regularization (LASSO/Elastic Net) on the voxel-specific amplitudes ($\beta_{k,v}$) associated with these numerous predictors. This reflects the neuroscientific assumption that, for any given voxel, only a subset of a large pool of potential continuous features (e.g., acoustic features, semantic embeddings, physiological recordings) will significantly drive its activity.

This extension aims to make `hrfals` a powerful tool for "feature-rich" fMRI modeling where data-driven HRF estimation is desired, leveraging and extending existing `hrfals` functionality while focusing on the "many continuous predictors" aspect and ensuring scalability and API cohesion.

## 2. Existing Model Framework and Proposed Extension

The underlying generative model remains:
$$ y_v(t) = \sum_{k=1}^{K} \beta_{k,v} (x_k \star h_v)(t) + \text{baseline}_v(t) + \varepsilon_v(t) $$

If $h$ is an FIR basis of length $d$, $(x_k \star h)(t)$ is the $t$-th element of $X_k h$, where $X_k$ is an $n \times d$ matrix of time-lagged versions of the $k$-th predictor $x_k(t)$. The model for voxel $v$ is:
$$ \mathbf{y}_v = \tilde{X}(h) \boldsymbol{\beta}_v + \text{baseline}_v + \boldsymbol{\varepsilon}_v $$
where $\tilde{X}(h) = [X_1 h, \dots, X_K h]$ is an $n \times K$ matrix of HRF-convolved predictors.

### Core Extension: Sparse $\beta$-Update

The estimation of $\boldsymbol{\beta}_v = [\beta_{1,v}, \dots, \beta_{K,v}]^\top$ for each voxel $v$ will be modified to solve an Elastic Net penalized regression:
$$ \hat{\boldsymbol{\beta}}_v = \arg\min_{\boldsymbol{\beta}_v} \frac{1}{2n} || \mathbf{y}_v^* - \tilde{X}(h) \boldsymbol{\beta}_v ||_2^2 + \lambda_{\beta,1} \left( (1-\alpha_\beta) \frac{1}{2} ||\boldsymbol{\beta}_v||_2^2 + \alpha_\beta ||\boldsymbol{\beta}_v||_1 \right) $$

where $\mathbf{y}_v^*$ is the BOLD signal after projecting out confounds and baseline model components. This will be implemented using `glmnet::glmnet` initially, with potential for a custom C++ coordinate descent solver for future optimization.

## 3. Key Algorithmic and Implementation Details

### 3.1. Predictor Standardization

- **Default: ON** via `standardize_predictors = TRUE` parameter
- Each original continuous predictor $x_k(t)$ will be z-scored (mean $\mu_k$, std $\sigma_k$)
- The $n \times d$ lagged matrix $X_k$ will be formed such that column $j$ (representing $x_k(t-(j-1)\Delta t)$) inherits this same scaling: $(x_k(t-(j-1)\Delta t) - \mu_k) / \sigma_k$
- This ensures all $d$ columns within $X_k$ (and thus all columns of $\tilde{X}(h)$) are on a comparable scale
- The $K$ means ($\mu_k$) and $K$ standard deviations ($\sigma_k$) will be stored in the returned `hrfals_fit` object
- Estimated $\hat{\beta}_{k,v}$ will be rescaled back to their original units before being returned

### 3.2. Formation of $X_k$ and Design Block Caching

- Input predictors $x_k(t)$ will be transformed into $n \times d$ lagged matrices $X_k$, where column $j$ represents $x_k(t-(j-1)\Delta t)$
- **Caching Control**: Parameter `cache_design_blocks = TRUE` controls whether the $K$ matrices $X_k$ are pre-computed and cached
- **Memory Management**: If caching is enabled but estimated memory footprint (approx. $K \times n \times d \times 8$ bytes) exceeds a practical threshold (e.g., 70% of available memory or 30-40 GB), the implementation will automatically switch to on-the-fly generation with a user message
- The $n \times K$ matrix $\tilde{X}(h)$ will be recomputed once per outer ALS iteration using the current $h$ estimate

### 3.3. Warm-Start for $\beta$-Step

- **Default: ON** via `beta_penalty$warm_start = TRUE`
- **Decision Rule**: 
  ```R
  perform_warm_start <- beta_penalty$warm_start && (beta_penalty$l1 > 0)
  ```
- **Implementation**:
  - In the first outer ALS iteration, the first $\beta$-update uses an unpenalized (or minimally L2-penalized) ridge solution as starting point
  - For subsequent iterations, $\beta$s from the previous iteration serve as the warm start
  - This improves convergence for sparse solvers

### 3.4. Penalty Parameter Integration

To maintain backward compatibility and provide flexible control:

**If `beta_penalty$l1 > 0` (Sparse/Elastic Net path is active):**
- The `glmnet` solver will be called with its main regularization strength parameter $\lambda$ set to `beta_penalty$l1` and mixing parameter set to `beta_penalty$alpha`
- The existing `lam_beta` argument provides an *additional* L2 penalty added to the Elastic Net's L2 component
- Effective L2 penalty becomes: `lambda_glmnet_L2_component + lam_beta`

**If `beta_penalty$l1 == 0` (Pure Ridge path):**
- The `beta_penalty$alpha` is ignored
- The L2 penalty strength is simply `lam_beta` (recovering current `hrfals` behavior)

### 3.5. $h$-Step (HRF Update) Scalability

The existing block-wise $h$-update logic in `cf_als_engine.R` will be largely preserved with optimizations for large $K$:

- **Diagonal Approximation**: When `fullXtX_flag = FALSE`, use diagonal approximation crucial for large $K$
- **Pre-computation**: The $K \times d$ matrix $s_{k,j} = \sum_t x_k(t-j\Delta t)^2$ will be pre-computed once
- **Efficient RHS Computation**: The RHS term $b_j = \sum_{k,v}\beta_{kv}\sum_{t}x_k(t-j\Delta t)r_{v}^*(t)$ will be computed to ensure linear scaling with $K$
- **Spatial Regularization**: Logic for handling spatial Laplacian (`lambda_s > 0`) and switching between direct solve and CG remains unchanged

## 4. API Modifications

The main user-facing wrapper `hrfals()` will be updated with grouped control parameters:

```R
hrfals(
    fmri_data_obj,
    event_model,
    hrf_basis,
    confound_obj = NULL,
    baseline_model = NULL,
    
    # Existing lambda arguments
    lam_beta = 10,               # L2 penalty strength
    lam_h = 1,                   # For HRF coefficients
    lambda_s = 0,                # For spatial smoothing of HRF
    
    # New grouped arguments for sparse beta estimation
    beta_penalty = list(
        l1 = 0,                  # L1 strength. >0 enables sparse betas
        alpha = 1,               # Elastic Net mix (1=LASSO, 0=Ridge)
        warm_start = TRUE        # Use ridge warm-start for sparse beta-step
    ),
    
    # New grouped arguments for design control
    design_control = list(
        standardize_predictors = TRUE,  # Z-score continuous x_k(t) predictors
        cache_design_blocks = TRUE      # Cache n x d lagged predictor matrices X_k
    ),
    
    # Other existing control parameters
    laplacian_obj = NULL,
    h_solver = c("direct", "cg", "auto"),
    cg_max_iter = 100,
    cg_tol = 1e-4,
    R_mat = NULL,
    fullXtX = FALSE,
    max_alt = 1,
    ...
)
```

**Convenience Alias**: `hrfals_sparse(...)` will be provided that calls `hrfals(...)` with default `beta_penalty$l1 > 0`.

## 5. Hyperparameter Tuning for Sparse Penalties

### Recommended Approach
- Global grid-search on a subset of voxels with limited CF-ALS iterations for selecting `beta_penalty$l1` and `beta_penalty$alpha`
- Default `alpha=1` (LASSO) is a good starting point; `alpha=0.5` can be explored if multicollinearity is high

### Helper Function
A tuning helper function `tune_beta_l1_hrfals()` will be provided:

```R
tune_beta_l1_hrfals(
    fmri_data_obj_subset,           # Y from a subset of voxels
    event_model,
    hrf_basis,
    l1_grid = 10^seq(-3, 0, length.out = 7),
    alpha_value = 1,                # Alpha to use during tuning
    n_outer_iterations_cfals = 1,   # For CF-ALS method
    other_hrfals_args = list(lam_beta=0.01, lam_h=0.01),
    cv_voxel_subset_train_prop = 0.7,
    seed = NULL
)
```

This function iterates through `l1_grid`, runs simplified `hrfals` on training data, predicts on test data, calculates MSE, and returns optimal hyperparameters.

## 6. Output (`hrfals_fit` object)

The `hrfals_fit` object will be enhanced to include:

- `h_coeffs` ($d \times v$) - HRF coefficients
- `beta_amps` ($K \times v$) - potentially sparse amplitudes, rescaled to original predictor units if standardized
- **Enhanced metadata**:
  - `beta_penalty` settings used for the fit
  - `predictor_means` and `predictor_sds` if `standardize_predictors = TRUE`
  - Expanded `lambdas` slot including sparse penalty parameters
- Existing fields: `method_used`, `call`, `fmrireg_hrf_basis_used`, `phi_recon_matrix`, `design_info`, `residuals`, `reconstructed_hrfs`, `gof_per_voxel`

## 7. Testing and Validation

### Unit Tests

**Backward Compatibility Test**: 
- Run existing `hrfals` test suite with `beta_penalty$l1 = 0` 
- Confirm numerically identical results to baseline runs

**Sparsity Recovery Test**:
- Simulate data with $K \approx 50-100$ predictors, where only ~5 have non-zero true $\beta$s
- Verify >90-95% of true zero $\beta_{k,v}$ are estimated as exactly zero
- Verify non-zero $\beta_{k,v}$ are recovered with reasonable accuracy

**Stress Test**:
- Test with $K=500, n=500, d=20$ for both `cache_design_blocks = TRUE/FALSE`
- Monitor RAM usage and runtime (target <3 min on 16-core machine)

**Feature Validation Tests**:
- Verify `standardize_predictors` and `warm_start` options function correctly
- Test graceful fallback when memory limits are exceeded

## 8. Documentation and Examples

### Comprehensive Documentation
- Update `hrfals` documentation to detail new parameters and their interactions
- Clearly document predictor standardization process and $\hat{\beta}_{k,v}$ reporting
- Document penalty parameter interactions and backward compatibility

### Demonstration Vignette
A vignette titled "Estimating Shared HRFs with Many Continuous Features (e.g., Movie fMRI)" will showcase:

1. Loading/simulating BOLD data and matrix of $K$ continuous predictors (semantic PCs, acoustic features)
2. Setting up `event_model` to incorporate these predictors
3. Running `hrfals()` with appropriate `beta_penalty` settings (using tuning helper)
4. Visualizing estimated shared HRF
5. Visualizing sparse $\beta_{k,v}$ matrix (heatmaps, profiles for selected voxels)
6. Interpreting which features drive which voxels

## 9. Conclusion

This extension significantly enhances the `hrfals` package by enabling the estimation of a shared HRF in the presence of many continuous predictors. By integrating sparse estimation for predictor amplitudes within the established CF-ALS framework, it offers a powerful and scalable solution for feature-rich fMRI modeling, while maintaining compatibility with `fmrireg`'s design specification tools. 

The core extension involves replacing the ridge $\beta$-step with an Elastic Net $\beta$-step and adding infrastructure for predictor standardization and design block caching, without requiring fundamental rewrites of existing `hrfals` engines. The pragmatic approach to caching, hyperparameter tuning, and graceful memory management should ensure both usability and computational feasibility for real-world neuroscience applications with hundreds of continuous predictors.

Okay, here's a granular, ticketed list of items to implement the "Sparse Beta Estimation for Shared HRF with Many Continuous Predictors" proposal within the `hrfals` package, aiming for accuracy, robustness, and minimal bugs.

This list assumes the `hrfals` package structure you provided (with `cf_als_engine.R`, `estimate_hrf_cfals.R`, etc.) is the starting point.

---

## Granular Implementation Plan: Sparse Beta CF-ALS for `hrfals`

**Phase 1: Core Engine Modifications & API Updates (`cf_als_engine.R`, `estimate_hrf_cfals.R`, Wrappers)**

*   **Ticket SP-CFALS-001: Verify/Ensure Block $h$-Update in `cf_als_engine`**
    *   **Task:** Double-check the existing $h$-update logic in `cf_als_engine.R`. Confirm it performs a block update for the entire $h$ vector (or $h_v$ per voxel if `lambda_s = 0`) simultaneously, not a coordinate-wise update for individual $h_j$ coefficients.
    *   **Acceptance:** Code review confirms block update. The RHS for the $h$-solve correctly uses the full relevant residual.
    *   **Note:** Based on prior review, this seems to be the case. This ticket is primarily for formal verification.

*   **Ticket SP-CFALS-002: Add New Control Parameters to `hrfals()` and `estimate_hrf_cfals()` Wrappers**
    *   **Task:** Modify `hrfals()` and the underlying `estimate_hrf_cfals()` (or whichever is the main internal dispatcher) to accept:
        *   `beta_penalty = list(l1 = 0, alpha = 1, warm_start = TRUE)`
        *   `design_control = list(standardize_predictors = TRUE, cache_design_blocks = TRUE)`
    *   **Acceptance:** Functions accept these new list arguments with specified defaults. Parameters are passed down to `cf_als_engine` or used in pre-processing steps.
    *   **Consideration:** Decide if existing `lam_beta` in `hrfals()` maps to an L2 component within `beta_penalty` or is a separate additive ridge. Proposal favors the latter for clarity (see SP-CFALS-006).

*   **Ticket SP-CFALS-003: Implement Predictor Standardization (`design_control$standardize_predictors`)**
    *   **Task:** In `estimate_hrf_cfals` (or its design preparation stage like `create_cfals_design`), if `standardize_predictors = TRUE`:
        1.  Identify continuous predictors $x_k(t)$ that are *not* already standardized basis functions (e.g., output of `fmrireg::Scale()`).
        2.  For each such $x_k(t)$, compute its mean ($\mu_k$) and standard deviation ($\sigma_k$) across time, ignoring NaNs. Handle zero/NA SD by using a small epsilon.
        3.  Standardize $x_k(t)' = (x_k(t) - \mu_k) / \sigma_k$.
        4.  Store $\mu_k$ and $\sigma_k$ (K-length vectors) in the `hrfals_fit` object (e.g., in `design_info` or a new slot like `predictor_scaling_info`).
        5.  Use $x_k(t)'$ for forming the lagged design blocks $X_k$.
    *   **Acceptance:** Standardization is applied correctly. Scaling factors are stored. $\hat{\beta}_{k,v}$ estimated on standardized predictors are correctly rescaled back to original units before being stored in `hrfals_fit$beta_amps`.

*   **Ticket SP-CFALS-004: Implement Caching for Lagged Design Blocks $X_k$ (`design_control$cache_design_blocks`)**
    *   **Task:** In the design preparation stage (e.g., `create_cfals_design` or early in `cf_als_engine`):
        1.  Estimate memory for all $K$ matrices $X_k$ (each $n \times d$).
        2.  If `cache_design_blocks = TRUE` AND estimated memory < threshold (e.g., 32GB or 70% of `memory.limit()`):
            *   Pre-compute and store all $X_k$ matrices.
        3.  Else (cache is FALSE or memory limit exceeded):
            *   Set an internal flag to generate $X_k$ or $X_k h$ components on-the-fly.
            *   If auto-disabled due to memory, issue a `message()`.
    *   **Acceptance:** Caching logic implemented. $X_k$ are either pre-computed or flag is set for on-the-fly generation.

*   **Ticket SP-CFALS-005: Implement Sparse $\beta$-Step in `cf_als_engine`**
    *   **Task:** Modify the $\beta$-update loop within `cf_als_engine`:
        1.  At the start of each $\beta$-step (once per outer ALS iteration), compute $\tilde{X}(h) = [X_1 h, \dots, X_K h]$ using current $h$ and (cached or on-the-fly) $X_k$.
        2.  **Warm-start logic:**
            *   If `beta_penalty$warm_start = TRUE` AND `beta_penalty$l1 > 0` AND it's the first outer ALS iteration:
                *   Perform an initial $\beta_v$ update for each voxel using ridge regression (e.g., `lambda_l1 = 0`, small `lam_beta` or existing L2 from `beta_penalty$alpha` if `l1` was temporarily zeroed). This result becomes the `beta_start_values_for_sparse_solver`.
            *   Else, `beta_start_values_for_sparse_solver` are the $\beta_v$ from the previous outer iteration.
        3.  For each voxel $v$:
            *   If `beta_penalty$l1 > 0`:
                *   Solve the Elastic Net problem for $\boldsymbol{\beta}_v$ using $\tilde{X}(h)$ and $y_v^*$. Use `glmnet::glmnet` with `lambda = beta_penalty$l1` (appropriately scaled if `glmnet` uses different convention), `alpha = beta_penalty$alpha`. Pass `beta_start_values_for_sparse_solver[,v]` if available and supported by solver.
            *   Else (`beta_penalty$l1 == 0`):
                *   Use the existing ridge regression $\beta$-update logic (currently in `cf_als_engine`), using `lam_beta` (from main args) as the ridge strength.
    *   **Acceptance:** $\beta$-step correctly switches between ridge and Elastic Net based on `beta_penalty$l1`. Warm-start is implemented. `glmnet` is called with correct parameters.

*   **Ticket SP-CFALS-006: Clarify and Implement Interaction of `lam_beta` and `beta_penalty`**
    *   **Task:** Define the precise calculation of L2 penalty in the $\beta$-step:
        *   If `beta_penalty$l1 == 0`: L2 penalty strength is `lam_beta`.
        *   If `beta_penalty$l1 > 0`: The Elastic Net solver (`glmnet`) handles its own L2 component via `alpha`. The `lam_beta` argument to `hrfals` will act as an *additional* ridge penalty summed with the Elastic Net's L2 term.
            *   `glmnet`'s effective $\lambda_{L2}$ is $\lambda \cdot (1-\alpha)/2$.
            *   The total L2 penalty on $\beta_v$ would be $\lambda_{\text{glmnet_L2_component}} + \text{lam_beta}$. This needs careful implementation if `glmnet` doesn't directly allow adding an external L2 penalty; might involve modifying the target function or applying two stages of regularization if simple addition isn't possible (less ideal).
            *   **Simpler alternative (closer to reviewer):** If `beta_penalty$l1 > 0`, `glmnet` is called with `lambda = beta_penalty$l1` and `alpha = beta_penalty$alpha`. The `lam_beta` argument is *ignored* by `glmnet` itself. If a user *really* wants an additional overall L2 shrinkage on top of what ElasticNet provides, this would be a feature for much later. For now, if `l1>0`, `beta_penalty$alpha` controls the L2/L1 mix for that primary penalty. If `l1=0`, then `lam_beta` controls a pure L2 penalty. This is the cleanest separation.
    *   **Acceptance:** The L2 penalty applied in the $\beta$-step is correctly determined and documented based on `lam_beta`, `beta_penalty$l1`, and `beta_penalty$alpha`. The simpler alternative is preferred.

*   **Ticket SP-CFALS-007: Ensure $h$-Step Scalability with Large $K$**
    *   **Task:** In `cf_als_engine`, when `fullXtX_flag = FALSE`:
        1.  Pre-compute the $K \times d$ matrix of scalars $s_{k,j} = \sum_t x_k(t-(j-1)\Delta t)^2$ using the (standardized, cached/on-the-fly) $X_k$ blocks.
        2.  Ensure the loop forming the diagonal Hessian elements $A_{v,jj}$ for $h_v$ (or $A_{jj}$ for global $h$) uses these $s_{k,j}$ and current $\beta_{kv}^2$ to maintain $O(K \cdot d)$ cost per voxel for this part.
        3.  Ensure the loop forming the RHS $b_{v,j}$ (or $b_j$) for $h_v$ (or global $h$) as $\sum_k \beta_{kv} (X_k^T r_v^* )_j$ also scales linearly with $K$. The pseudocode `for (k in 1:K) { xkj <- X_shifted[[k]][, j]; rhs_j += crossprod(xkj, resid_mat %*% B[k, ]) }` needs careful implementation of `resid_mat` and `B[k,]` to achieve this. `resid_mat` is $n \times v$, `B[k,]` is $1 \times v$. `resid_mat %*% B[k,]` is $n \times 1$.
    *   **Acceptance:** $h$-update computations involving sums over $K$ predictors are confirmed to be $O(K)$ and not $O(K^2)$.

**Phase 2: Output, Documentation, and Higher-Level Wrappers**

*   **Ticket SP-CFALS-008: Update `hrfals_fit` Object Output**
    *   **Task:** Modify the `hrfals_fit()` constructor and the returned object:
        1.  Store the full `beta_penalty` list used in a slot like `lambdas$beta_penalty`.
        2.  If `standardize_predictors = TRUE`, store the K-length vectors `predictor_means` and `predictor_sds` (e.g., in `design_info` or a new slot like `predictor_scaling_info`).
    *   **Acceptance:** `hrfals_fit` object contains the new penalty information and scaling factors.

*   **Ticket SP-CFALS-009: Create Alias `hrfals_sparse()`**
    *   **Task:** Implement the thin wrapper:
        ```R
        hrfals_sparse <- function(..., beta_penalty = list(l1 = 0.05, alpha = 1, warm_start = TRUE), design_control = list(standardize_predictors = TRUE, cache_design_blocks = TRUE)) {
            hrfals(..., beta_penalty = beta_penalty, design_control = design_control)
        }
        ```
        Export this function.
    *   **Acceptance:** Alias function created and works as expected.

*   **Ticket SP-CFALS-010: Update Documentation**
    *   **Task:** Thoroughly document the new parameters in `hrfals()` (`beta_penalty`, `design_control`) and their sub-options.
        *   Explain the interaction between `lam_beta` and `beta_penalty`.
        *   Explain predictor standardization and how reported $\beta$s are in original units.
        *   Explain the `cache_design_blocks` mechanism and its implications.
        *   Describe the warm-start behavior.
        *   Provide guidance on choosing `beta_penalty$l1` and `beta_penalty$alpha`.
    *   **Acceptance:** All new functionalities and parameters are clearly documented. Man pages updated.

*   **Ticket SP-CFALS-011: Add NEWS.md Entry**
    *   **Task:** Add a bullet point under "New Features" in `NEWS.md` announcing "Sparse CF-ALS: Estimation of a shared HRF with hundreds of continuous predictors via LASSO/Elastic Net regularization on beta amplitudes."
    *   **Acceptance:** NEWS.md updated.

**Phase 3: Advanced Features & Testing**

*   **Ticket SP-CFALS-012: Develop Helper `tune_beta_l1_hrfals()` (Optional - Post MVP)**
    *   **Task:** Implement the hyperparameter tuning helper function as sketched in prior discussions, allowing users to find a data-driven `beta_penalty$l1`.
    *   **Acceptance:** Helper function implemented, tested, and documented.

*   **Ticket SP-CFALS-013: Implement Unit Tests - Core Functionality**
    *   **Task:**
        1.  Test that `hrfals` with `beta_penalty$l1 = 0` reproduces results from the current `hrfals` (before this extension) to a high tolerance (regression test).
        2.  Create a simple simulation with $K$ predictors, where only a few have true non-zero $\beta$s. Test that with an appropriate `beta_penalty$l1 > 0`, the sparse $\beta$-step correctly identifies these and shrinks others to zero.
        3.  Test `standardize_predictors = TRUE` vs. `FALSE` and verify $\beta$ rescaling.
        4.  Test `beta_penalty$warm_start = TRUE` vs. `FALSE` and check for convergence differences (though exact equivalence isn't expected, both should converge to similar sparse solutions).
    *   **Acceptance:** All core functionality tests pass.

*   **Ticket SP-CFALS-014: Implement Unit Test - Stress Test**
    *   **Task:** Create a unit test with $K=500$ predictors, $n=500$ timepoints, $d=20$ FIR lags.
        1.  Run with `design_control$cache_design_blocks = TRUE`. Monitor memory (manually or with tools if CI supports).
        2.  Run with `design_control$cache_design_blocks = FALSE`. Measure runtime and assert it's within a reasonable target (e.g., <3 minutes on a reference CI machine or typical developer laptop).
    *   **Acceptance:** Stress tests complete within performance/memory targets.

*   **Ticket SP-CFALS-015: Develop Vignette for "Many Continuous Predictors"**
    *   **Task:** Create a new vignette demonstrating the use case:
        1.  Simulating/loading fMRI data and a matrix of $K \approx 100-300$ continuous features.
        2.  Setting up the `fmrireg::event_model` to use these features (e.g., perhaps each $x_k(t)$ is a column in the `data` frame passed to `event_model`, and the formula is `onset ~ hrf(x1) + hrf(x2) + ...`, or a more programmatic way to specify many `hrf(x_k)` terms if `fmrireg` supports it, or if this CF-ALS variant will take `X_list_raw` directly). *This part needs careful thought on how users will specify many $x_k(t)$ to `hrfals`.*
        3.  Running `hrfals()` with the new sparse beta estimation.
        4.  Visualizing the estimated shared HRF and the sparse $\beta_{k,v}$ weights (e.g., heatmap of $\beta_v$, or which features are selected for example voxels).
    *   **Acceptance:** Vignette is clear, runnable, and effectively showcases the new functionality.

**Considerations for `fmrireg::event_model` interaction (re: Vignette task SP-CFALS-015):**

The current `fmrireg::event_model` using formulas like `onset ~ hrf(x1) + hrf(x2) + ... + hrf(xK)` would create $K$ separate `event_term` objects. The current `estimate_hrf_cfals` targets *one* `event_term`.

*   **Option A:** The user provides a *single* `event_term` to `estimate_hrf_cfals` where this term is constructed using `fmrireg::Ident(x1, x2, ..., xK)`. Then `fmrireg::create_fmri_design` (which `hrfals` calls via `create_cfals_design`) would need to correctly produce $K$ separate $X_k$ design blocks from this single `Ident` term when combined with the `hrf_basis_for_cfals`. This is the most `fmrireg`-idiomatic way.
*   **Option B:** `estimate_hrf_cfals` is modified to accept a *list* of $K$ continuous predictor time series $x_k(t)$ directly, bypassing `event_model` for this specific use case, and forms the $X_k$ internally. This is less integrated but might be simpler if `fmrireg::Ident()` doesn't naturally produce the required $X_k$ structure for CF-ALS.

Option A is preferable for consistency. The `create_fmri_design` function (from `hrfals`) would need to be smart about how it processes an `event_model` term created with `hrf(Ident(x1, ..., xK), basis = hrf_basis_for_cfals)` to yield the list of $K$ separate $n \times d$ matrices. This seems achievable as `fmrireg::convolve.event_term` would produce an $n \times (K \cdot d)$ matrix, which could then be split.

This ticketed list should provide a solid, actionable plan.