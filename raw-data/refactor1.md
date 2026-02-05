

**Revised Proposal: Robust & Flexible CF-ALS Suite for `hrfals` (v3.1 - `fmrireg` Consumption Detailed)**

**1. Objective:** (Same as "All-R Refactor v2.1")
To develop an all-R CF-ALS suite in the `hrfals` package, consuming `fmrireg` objects, for data-driven HRF estimation, offering `ls_svd_only`, `ls_svd_1als` (default), and iterative `cf_als` modes.

**2. Core Algorithms & Modules (All R Implementation within `hrfals` package):**

*   **A. CF-ALS Design & Projection Utilities (in `hrfals/R/cfals_design_utils_hrfals.R`):**
    *   **Function:** `hrfals::prepare_cfals_inputs_from_fmrireg_term(fmri_data_obj, fmrireg_event_model, target_event_term_name, fmrireg_hrf_basis, confound_obj, hrf_shape_duration_sec, hrf_shape_sample_res_sec)`
        *   **Inputs:**
            *   `fmri_data_obj`: `fmridesign::fmri_dataset` or BOLD matrix.
            *   `fmrireg_event_model`: An `fmridesign::event_model` object.
            *   `target_event_term_name`: Character string, the name of the *single* `event_term` within `fmrireg_event_model` to process with CF-ALS.
            *   `fmrireg_hrf_basis`: An `fmridesign::HRF` object to be used for *this CF-ALS fit*. **Must have `fmridesign::nbasis(fmrireg_hrf_basis) > 1`.**
            *   `confound_obj`: Optional confound matrix.
            *   `hrf_shape_duration_sec`, `hrf_shape_sample_res_sec`: For `Phi` and `h_ref_shape`.
        *   **Steps:**
            1.  Extract `Y_raw`, `Z_confounds`. Handle NA rows in `Y_raw`.
            2.  **Select Target `event_term` from `fmrireg_event_model`:**
                *   `fmrireg_directive:`
                    ```R
                    all_terms <- fmridesign::terms(fmrireg_event_model)
                    if (!target_event_term_name %in% names(all_terms)) {
                      stop("Target event term '", target_event_term_name, "' not found in event_model.")
                    }
                    target_term_obj <- all_terms[[target_event_term_name]]
                    if (!inherits(target_term_obj, "event_term")) {
                        stop("Selected term is not an 'event_term', CF-ALS requires an event_term.")
                    }
                    ```
            3.  **Get Unconvolved Predictors for the Target Term (`X_unconvolved_term` - an `n x k_term_conditions` matrix):**
                *   `fmrireg_directive:`
                    ```R
                    # design_matrix.event_term gives the unconvolved predictors (e.g., stick functions per condition level)
                    X_unconvolved_term <- fmridesign::design_matrix(target_term_obj, drop.empty = TRUE) 
                    # drop.empty = TRUE is important to match number of estimable conditions
                    ```
                *   `k_conditions = ncol(X_unconvolved_term)`. If `k_conditions == 0`, stop or warn.
            4.  `d_basis_dim = fmridesign::nbasis(fmrireg_hrf_basis)`.
                *   If `d_basis_dim <= 1`, `stop("CF-ALS requires an hrf_basis with nbasis > 1.")`.
            5.  **Generate `X_list_raw` (`k_conditions` list of `n x d_basis_dim` matrices):**
                *   For each condition `c` (column `c`) of `X_unconvolved_term`:
                    *   `raw_onsets_c_timeseries = X_unconvolved_term[, c]` (an `n x 1` vector representing that condition's unconvolved predictor).
                    *   `X_list_raw[[c]]` (`n x d_basis_dim`):
                        *   Column `j` is `convolve_timeseries_with_single_basis(raw_onsets_c_timeseries, fmrireg_hrf_basis, basis_function_index = j, sampling_frame = fmrireg_event_model$sampling_frame)`.
                        *   *(This helper `convolve_timeseries_with_single_basis` needs to be implemented in `hrfals`. It would take an `fmridesign::HRF` object, extract its `j`-th basis function (e.g., by creating a temporary `HRF` object for it), and perform convolution, possibly using `fmridesign::evaluate` and `stats::convolve` or a C++ helper.)*
                *   Apply NA row zeroing to `X_list_raw` consistent with `Y_raw`.
            6.  **Generate `Phi_recon_matrix` and `h_ref_shape_canonical_p_dim` (within `hrfals`):**
                *   `fmrireg_directive (for Phi):`
                    ```R
                    time_points_for_shape = seq(0, hrf_shape_duration_sec, by = hrf_shape_sample_res_sec)
                    # Phi_recon_matrix will be p x d_basis_dim
                    Phi_recon_matrix = fmridesign::evaluate(fmrireg_hrf_basis, time_points_for_shape)
                    if (is.vector(Phi_recon_matrix) && d_basis_dim == 1) Phi_recon_matrix <- matrix(Phi_recon_matrix, ncol=1)
                    if (ncol(Phi_recon_matrix) != d_basis_dim) stop("Mismatch in Phi_recon_matrix columns and nbasis.")
                    ```
                *   Generate `p`-dim Glover shape for `h_ref_shape_canonical_p_dim`.
            7.  Confound Projection: `Y_proj`, `X_list_proj`.
        *   **Outputs:** `list(Y_proj, X_list_proj, d_basis_dim, k_conditions, Phi_recon_matrix, h_ref_shape_canonical_p_dim, n_timepoints, v_voxels, bad_row_idx, target_term_condition_names = colnames(X_unconvolved_term))`.

*   **B. Refactored `ls_svd_engine` (R Function - in `hrfals`):** (As before, uses outputs of `prepare_cfals_inputs_from_fmrireg_term`).

*   **C. Refactored `cf_als_engine` (R Function - iterative core, in `hrfals`):**
    *   **Input `R_mat_user`:** If `NULL`, `R_mat_eff = diag(d_basis_dim)`. *`hrfals` will not attempt to call a `fmridesign::penalty_matrix` S3 method unless `hrfals` itself defines one for `fmridesign::HRF` objects.*
    *   (Rest as before).

*   **D. Refactored `ls_svd_1als_engine` (R Function - in `hrfals`):** (Calls refactored engines).

*   **E. Main User-Facing Wrapper `hrfals::estimate_hrf_cfals(...)` (R):**
    1.  **Inputs:**
        *   `fmri_data_obj`: `fmridesign::fmri_dataset` or BOLD matrix.
        *   `fmrireg_event_model`: An `fmridesign::event_model` object.
        *   `target_event_term_name`: Character string (name of `event_term` in `fmrireg_event_model`).
        *   `hrf_basis_for_cfals`: An `fmridesign::HRF` object (e.g., `fmridesign::HRF_BSPLINE()`), **must have `nbasis > 1`**.
        *   (Other params: `confound_obj`, `method`, lambdas, `R_mat` (user custom), `fullXtX`, `max_alt`, `hrf_shape_duration`, `hrf_shape_resolution`).
    2.  Calls `prepare_cfals_inputs_from_fmrireg_term`.
    3.  Dispatches to R engines.
    4.  Post-processing:
        *   Uses `target_term_condition_names` from `prepare_cfals_inputs_from_fmrireg_term` for `rownames(fit$beta)`.
        *   Corrected prediction for residuals/R².
    5.  Return `hrfals_fit` object.

*   **F. Output Class & Methods (in `hrfals`):**
    *   `hrfals_fit` class: Store `fmrireg_hrf_basis_used` (the one passed to `estimate_hrf_cfals`).

**4. Sprint Plan (All-R Refactor in `hrfals` Package, Focusing on Single `fmridesign::event_term`):**

**Sprint 0: `hrfals` Package Setup & Core `fmrireg` Interaction Utilities**
*   **`HRFALS-V2.1-S0-T1`: Package Structure & Dependencies** (Imports `fmrireg`).
*   **`HRFALS-V2.1-S0-T2`: Implement `hrfals:::convolve_timeseries_with_single_basis` Helper:**
    *   **Task:** Takes `raw_timeseries (n x 1)`, `fmrireg_hrf_obj`, `basis_function_index (j)`, `fmrireg_sampling_frame`. Creates a temporary single-basis `fmridesign::HRF` from `fmrireg_hrf_obj[[j]]` and convolves.
    *   **DoD:** Helper correctly convolves a timeseries with the j-th basis function of any `fmridesign::HRF` object.
*   **`HRFALS-V2.1-S0-T3`: Implement `hrfals::prepare_cfals_inputs_from_fmrireg_term` Utility:**
    *   Implement logic to select `target_term_obj` from `fmrireg_event_model`.
    *   Use `fmridesign::design_matrix(target_term_obj)` to get `X_unconvolved_term`.
    *   Loop `k_conditions` times, calling `convolve_timeseries_with_single_basis` `d_basis_dim` times to form `X_list_raw`.
    *   Implement NA row handling, `Phi_recon_matrix` generation (using `fmridesign::evaluate`), canonical `h_ref_shape_p_dim` generation, confound projection.
    *   **DoD:** Utility correctly prepares all inputs for R engines from a *single specified term* within an `fmridesign::event_model`. Unit tested.
*   **`HRFALS-V2.1-S0-T4`: Define `hrfals_fit` S3 Class & Basic Print/Summary.**

**Sprint 1: Refactor `ls_svd_engine` in `hrfals` (No change from previous sprint plan `CFALS-R2.1-S1-T1`, `CFALS-R2.1-S1-T2`)**

Ticket CFALS-R2.1-S2-T1: Refactor cf_als_engine - Precomputation & XtY Strategy
Current State (from repomix-output.xml R/cf_als_engine.R):
Precomputes XtX_list = lapply(X_list_proj, crossprod).
Precomputes XtY_list = lapply(X_list_proj, function(X) crossprod(X, Y_proj)). This creates k matrices, each d x v.
If fullXtX_flag, precomputes XtX_full_list (all X_lᵀX_m).
Task Details:
Modify cf_als_engine_r Signature:
Add a new argument precompute_xty_flag (boolean, default TRUE).
The inputs X_list_proj (k list of n x d matrices) and Y_proj (n x v matrix) are received from the calling wrapper (e.g., hrfals::estimate_hrf_cfals via prepare_cfals_inputs_from_fmrireg_term).
Implement XtY Strategy based on precompute_xty_flag:
If precompute_xty_flag = TRUE (Default):
Retain the current logic: XtY_list = lapply(X_list_proj, function(X) crossprod(X, Y_proj)).
Memory Warning (New): Calculate the approximate size of XtY_list (e.g., k_conditions * d_basis_dim * v_voxels * 8 bytes). If it exceeds a threshold (e.g., 1-2 GB, configurable or hardcoded with a message), print a message() suggesting precompute_xty_flag = FALSE for large datasets.
If precompute_xty_flag = FALSE:
Do not compute the full XtY_list.
The per-voxel loops for β-update and h-update (in subsequent tickets) will need to compute XtY_vx (a d x 1 vector for the current voxel vx for each condition c: crossprod(X_list_proj[[c]], Y_proj[,vx])) on the fly.
Precomputation of XtX terms:
Maintain precomputation of XtX_list (list of k d x d matrices: X_cᵀX_c).
Maintain precomputation of XtX_full_list (list of k x k d x d matrices: X_lᵀX_m) if fullXtX_flag = TRUE.
Inputs to this specific ticket's refactored part (within cf_als_engine_r): X_list_proj, Y_proj, d_basis_dim, k_conditions, v_voxels, fullXtX_flag, precompute_xty_flag.
Outputs from this specific ticket's refactored part: XtX_list, XtX_full_list (if applicable), and either full XtY_list or just the components needed for on-the-fly calculation (X_list_proj, Y_proj passed through).
Definition of Done:
cf_als_engine_r correctly implements both XtY precomputation strategies.
Memory warning for full XtY_list precomputation is implemented.
XtX terms are precomputed as before.
Unit tests verify both XtY strategies (e.g., by checking intermediate values are the same if computed on-the-fly vs. precomputed for a small case).
Ticket CFALS-R2.1-S2-T2: Refactor cf_als_engine - β and h Update Loops
Current State (from R/cf_als_engine.R):
Per-voxel R loops for β and h updates.
β-update forms G_vx (k x k) and solves.
h-update forms lhs_vx (d x d) and rhs_vx (d x 1), solves. lambda_h currently applies to an identity matrix if R_mat is NULL.
Task Details:
Modify cf_als_engine_r Signature:
Accept R_mat_eff (a d x d penalty matrix, e.g., diag(d_basis_dim) for ridge, or a custom smoothness matrix). This is passed from hrfals::estimate_hrf_cfals.
Refactor Per-Voxel β-update Loop:
Iterate vx from 1 to v_voxels.
h_vx = H_current[,vx].
XtY_vx Handling:
If precompute_xty_flag = TRUE: XtY_c_vx = XtY_list[[c]][,vx].
If precompute_xty_flag = FALSE: XtY_c_vx = crossprod(X_list_proj[[c]], Y_proj[,vx]).
DhTy_vx = vapply(seq_len(k_conditions), function(c) crossprod(h_vx, XtY_c_vx_for_condition_c), numeric(1)).
Form G_beta_matrix_vx (k x k) using h_vx, XtX_list, XtX_full_list (if fullXtX_flag), correctly as per reviewed C++ sketch (i.e., G_vx[l,m] = h_vxᵀ (X_lᵀX_m) h_vx).
Solve (G_beta_matrix_vx + lambda_b * diag(k_conditions)) β = DhTy_vx for B_current[,vx], using robust_chol_solve.
Refactor Per-Voxel h-update Loop:
Iterate vx from 1 to v_voxels.
b_vx = B_current[,vx].
XtY_vx Handling: (Same as in β-update).
rhs_h_vx = numeric(d_basis_dim). Loop l from 1:k_conditions, rhs_h_vx = rhs_h_vx + b_vx[l] * XtY_l_vx_for_condition_l.
lhs_h_matrix_vx = lambda_h * R_mat_eff. (Note: lambda_h is now a scalar multiplier for the entire R_mat_eff).
Loop l from 1:k_conditions (and m if fullXtX_flag):
lhs_h_matrix_vx = lhs_h_matrix_vx + b_vx[l]*b_vx[m]* (if fullXtX XtX_full_list[[l,m]] else if l==m XtX_list[[l]] else 0_matrix).
Solve lhs_h_matrix_vx h = rhs_h_vx for H_current[,vx], using robust_chol_solve.
Definition of Done:
β-update and h-update loops correctly implemented in R.
Both XtY strategies (precompute_xty_flag = TRUE/FALSE) are functional within the loops.
R_mat_eff is correctly incorporated into the h-update.
Unit tests verify the correctness of a single β and h update step for a small case, with and without fullXtX_flag, and with a non-identity R_mat_eff.
Ticket CFALS-R2.1-S2-T3: Refactor cf_als_engine - Convergence & Final Identifiability (on Shapes)
Current State (from R/cf_als_engine.R):
ALS loop runs for fixed max_alt iterations. No explicit convergence check within the loop.
Final identifiability operates on h_current (coefficients) using a d-dimensional h_ref_shape_norm.
Task Details:
Modify cf_als_engine_r Signature:
Accept Phi_recon_matrix (p x d) and h_ref_shape_canonical_p_dim (p x 1) from prepare_cfals_inputs_from_fmrireg_term.
Accept tol_als (convergence tolerance).
Implement Convergence Check within ALS Loop:
Inside the iter_als loop (after h-update):
Calculate max_abs_diff_B = max(abs(B_current - B_old)).
Calculate max_abs_diff_H = max(abs(H_current - H_old)).
If max_abs_diff_B < tol_als * (1 + max(abs(B_old))) AND max_abs_diff_H < tol_als * (1 + max(abs(H_old))), then break the loop.
Store iterations_run = iter_als.
Implement Final Identifiability (on Reconstructed Shapes):
After the ALS loop (or if max_alt = 0, applied to H_init, B_init):
H_to_constrain = H_current, B_to_constrain = B_current.
reconstructed_H_iter = Phi_recon_matrix %*% H_to_constrain (p x v).
scl_factors = apply(abs(reconstructed_H_iter), 2, max).
alignment_scores = colSums(reconstructed_H_iter * h_ref_shape_canonical_p_dim).
Determine flip_signs.
effective_scl = pmax(scl_factors, epsilon_scale).
H_final = sweep(H_to_constrain, MARGIN = 2, STATS = flip_signs / effective_scl, FUN = "*").
B_final = sweep(B_to_constrain, MARGIN = 2, STATS = flip_signs * effective_scl, FUN = "*").
Zero out H_final and B_final columns where scl_factors was negligible.
Return Values:
list(h_coeffs=H_final, beta_amps=B_final, iterations_run=iterations_run, reconstructed_hrfs=Phi_recon_matrix %*% H_final).
Definition of Done:
ALS loop includes a robust convergence check.
Final identifiability is performed on reconstructed HRF shapes using Phi_recon_matrix.
Engine returns all specified outputs.
Unit tests verify convergence behavior and identifiability with a non-FIR basis (where Phi_recon_matrix is not identity).
Ticket CFALS-R2.1-S2-T4: Refactor ls_svd_1als_engine_r
Current State (from R/ls_svd_1als_engine.R):
Calls old ls_svd_engine.
Re-implements one iteration of β and h updates (similar to cf_als_engine but for one pass).
Performs its own identifiability.
Task Details:
Modify ls_svd_1als_engine_r to be a thin wrapper:
Call the (now refactored) ls_svd_engine_r to get H_init, B_init, and Gamma_hat.
Call the (now refactored) cf_als_engine_r with H_init, B_init, relevant parameters (lambdas, flags, Phi_recon_matrix, h_ref_shape_canonical_p_dim), and crucially, max_als_iter = 1.
The final identifiability is now handled inside cf_als_engine_r.
Return Values:
list(h=fit_cf_als$h_coeffs, beta=fit_cf_als$beta_amps, reconstructed_hrfs=fit_cf_als$reconstructed_hrfs, h_ls_svd=H_init, beta_ls_svd=B_init, iterations_run_als=fit_cf_als$iterations_run).
Definition of Done:
ls_svd_1als_engine_r correctly uses the refactored ls_svd_engine_r and cf_als_engine_r.
Its output structure is consistent.
Passes tests ensuring it's equivalent to cf_als_engine_r with max_alt=1.

**Sprint 3: Top-Level Wrapper & Full `fmrireg` Integration Testing**
*   **`HRFALS-V2.1-S3-T1`: Implement `hrfals::estimate_hrf_cfals` User-Facing Wrapper:**
    *   Takes `fmrireg_event_model` and `target_event_term_name`.
    *   Takes `hrf_basis_for_cfals` (an `fmridesign::HRF` object, check `nbasis > 1`).
    *   Calls `prepare_cfals_inputs_from_fmrireg_term`.
    *   Handles `penalty_R_mat_type` to generate simple `R_mat` or use custom.
    *   Corrected prediction formula for residuals/R².
    *   **DoD:** Main user function operational for a single target term.
*   **`HRFALS-V2.1-S3-T2`: Enhanced `hrfals_fit` Output Object & Plot Method:**
    *   Store `fmrireg_hrf_basis_used`, `target_event_term_name`.
    *   `plot()` uses `Phi_recon_matrix` from `prepare_cfals_inputs...` results.
    *   **DoD:** Output object comprehensive.
*   **`HRFALS-V2.1-S3-T3`: Comprehensive Integration Testing with `fmrireg` Objects:**
    *   Test `estimate_hrf_cfals` with various `fmridesign::HRF` types for `hrf_basis_for_cfals` (e.g., `HRF_SPMG1` (if nbasis made >1 for testing), `HRF_SPMG3`, `HRF_BSPLINE`), targeting different terms.
    *   **DoD:** High confidence in robust integration.

**Sprint 4: Documentation, Examples & R Profiling (As before, for `hrfals` package)**

This revised plan makes the interaction with `fmrireg` much more explicit and correct, particularly around how the "unconvolved design per basis function" (`X_list_raw`) is obtained when CF-ALS is to operate on a specific term from a larger `fmrireg` model, using a separate `fmrireg`-defined HRF basis for the estimation.

addendum:

You're right to ask for precision here. This is a fundamental point for how CF-ALS (or any method that *estimates* HRF coefficients from data) interacts with the stimulus design when the HRF is represented by a basis set.

Let's break it down:

**Goal of CF-ALS (or similar HRF estimation methods):**

We are trying to fit a model like this for each voxel `v`:

`Y_v ≈ sum_{c=1 to k} (X_raw_onsets_c ⋆ h_v) * β_{c,v}`

Where:
*   `Y_v` is the `n x 1` BOLD time series for voxel `v`.
*   `X_raw_onsets_c` is an `n x 1` vector representing the raw onsets for condition `c` (e.g., a series of delta functions at event times, or boxcars for block designs, *before any HRF convolution*).
*   `h_v` is the `p x 1` HRF shape (sampled at `p` time points) for voxel `v`.
*   `⋆` denotes convolution.
*   `β_{c,v}` is the scalar amplitude for condition `c` in voxel `v`.

Now, if we represent the unknown HRF `h_v` using a linear basis set of `d` basis functions `Φ = [φ_1, φ_2, ..., φ_d]` (where `Φ` is a `p x d` matrix, and each `φ_j` is a `p x 1` basis shape vector), then:

`h_v = Φ @ coeffs_h_v`

where `coeffs_h_v` is a `d x 1` vector of coefficients for voxel `v` that we want to estimate.

Substituting this into the model:

`Y_v ≈ sum_{c=1 to k} (X_raw_onsets_c ⋆ (Φ @ coeffs_h_v)) * β_{c,v}`

Due to the linearity of convolution:
`(A ⋆ (B @ C)) = (A ⋆ sum_j B_j * C_j) = sum_j (A ⋆ B_j) * C_j`
So, `(X_raw_onsets_c ⋆ (Φ @ coeffs_h_v)) = sum_{j=1 to d} (X_raw_onsets_c ⋆ φ_j) * coeffs_h_v[j]`

Let `X_convolved_basis_cj = (X_raw_onsets_c ⋆ φ_j)`. This is an `n x 1` time series: the raw onsets of condition `c` convolved with the `j`-th HRF basis function.

Then the sum becomes: `[X_convolved_basis_c1, X_convolved_basis_c2, ..., X_convolved_basis_cd] @ coeffs_h_v`.
Let `X_list_proj[[c]]` be this `n x d` matrix: `X_list_proj[[c]] = [X_convolved_basis_c1, ..., X_convolved_basis_cd]`.

So the model is:
`Y_v ≈ sum_{c=1 to k} (X_list_proj[[c]] @ coeffs_h_v) * β_{c,v}`

This is exactly the form `Y_v ≈ sum_c (X_c' @ h) * β_c` that CF-ALS uses, where:
*   `X_c'` is our `X_list_proj[[c]]`.
*   `h` is our `coeffs_h_v`.
*   `β_c` is our `β_{c,v}`.

**The "Unconvolved Design Per Basis Function" - Precise Form:**

What I mean by this refers to the matrix `X_list_proj[[c]]` (or its raw, pre-confound-projection version, `X_list_raw[[c]]`).

For a given **condition `c`** and an **HRF basis set `Φ` consisting of `d` basis functions `[φ_1, φ_2, ..., φ_d]`**:

The matrix `X_list_raw[[c]]` is an `n x d` matrix where:

*   **`n`** is the number of time points (scans) in the run/session.
*   **`d`** is the number of basis functions in the chosen HRF basis set.

*   **The `j`-th column of `X_list_raw[[c]]` is the time series resulting from convolving the raw stimulus onset function for condition `c` with the `j`-th HRF basis function `φ_j`.**

**Let's be more concrete with an example:**

Suppose:
*   `TR = 2s`.
*   `n = 100` scans (total time 200s).
*   Condition `A` has onsets at `t = 10s, 50s, 90s`.
*   We choose an HRF basis set with `d=3` basis functions: `φ_1(t)`, `φ_2(t)`, `φ_3(t)`.
    *   `φ_1(t)` might be a canonical Glover HRF.
    *   `φ_2(t)` might be its temporal derivative.
    *   `φ_3(t)` might be its dispersion derivative.
    (Each `φ_j(t)` is itself a shape defined over some duration, say `p=16` TRs = 32s).

To construct `X_list_raw[["A"]]` (the `100 x 3` matrix for Condition A):

1.  **Create the raw onset vector for Condition A:**
    `raw_onsets_A` = an `n x 1` vector (or a higher-resolution representation if event timing is not TR-locked). It would be 1 at scan indices corresponding to 10s, 50s, 90s, and 0 otherwise. (If events have duration, this would be a boxcar).

2.  **Column 1 of `X_list_raw[["A"]]`:**
    `X_list_raw[["A"]][, 1] = convolve(raw_onsets_A, φ_1)`
    This column represents the BOLD signal predicted by all onsets of condition A if the HRF was exactly `φ_1`.

3.  **Column 2 of `X_list_raw[["A"]]`:**
    `X_list_raw[["A"]][, 2] = convolve(raw_onsets_A, φ_2)`
    This column represents the BOLD signal predicted by all onsets of condition A if the HRF was exactly `φ_2`.

4.  **Column 3 of `X_list_raw[["A"]]`:**
    `X_list_raw[["A"]][, 3] = convolve(raw_onsets_A, φ_3)`

This `X_list_raw[["A"]]` is then what CF-ALS uses (after confound projection to become `X_list_proj[["A"]]`). When CF-ALS estimates `coeffs_h_v` for voxel `v`, these coefficients are the weights for these three pre-convolved basis function regressors. The estimated HRF shape for that voxel is then `h_v = coeffs_h_v[1]*φ_1 + coeffs_h_v[2]*φ_2 + coeffs_h_v[3]*φ_3`.

**Why this is different from a standard GLM design matrix with ONE HRF:**

*   In a standard GLM where you *assume* a fixed canonical HRF (say, `φ_canonical`), the design matrix for condition `c` would be a single `n x 1` vector: `X_canonical_c = convolve(raw_onsets_c, φ_canonical)`. The GLM estimates one beta `β_c` for this regressor.
*   In CF-ALS (and other HRF estimation methods using a basis set), we are trying to *estimate* the HRF shape. We do this by creating `d` regressors for *each condition*, where each regressor is that condition's onsets convolved with *one of the basis functions*. The GLM part of LS+SVD then estimates `d` coefficients *for each condition*. These `k` sets of `d` coefficients form the `G_v` matrix.

**How `fmrireg` typically handles this:**

The `fmrireg` package, when you specify `hrf(condition, basis = my_hrf_basis_obj)`, and `my_hrf_basis_obj` has `d` basis functions, `fmridesign::design_matrix(event_model)` usually returns a final design matrix where, for each level of `condition`, there are `d` columns. Each of these `d` columns is the result of convolving that condition level's onsets with one of the `d` basis functions from `my_hrf_basis_obj`.

If `fmridesign::design_matrix(event_model)` returns an `n x (k*d_total_effective)` matrix where `d_total_effective` might be `d_hrf_basis * d_other_parametric_modulators`, then:
*   `Xbig` in CF-ALS is this exact matrix (after confound projection).
*   `Gamma_hat` are the coefficients for these `k*d_total_effective` regressors.
*   Reshaping `Gamma_hat` to `G_v` (`d_total_effective x k`) correctly aligns these for the SVD if `d_total_effective` represents the "shape" part and `k` is the "amplitude" part for each condition.

The critical point is that `X_list_proj[[c]]` for CF-ALS *must* be the `n x d` matrix formed by `[ (onsets_c ⋆ φ_1), (onsets_c ⋆ φ_2), ..., (onsets_c ⋆ φ_d) ]`. If `fmrireg`'s `design_matrix` for a term involving `hrf(condition_factor, basis=hrf_object_d_dim)` produces columns like `conditionA_basis1, conditionA_basis2, ..., conditionB_basis1, ...`, then `Xbig` is already correctly formed by `cbind`-ing these terms. The `X_list_proj` is just a conceptual grouping of these columns per original condition.

So, "unconvolved design per basis function" means "design matrix where each column is (stimulus onsets for condition `c`) convolved with (HRF basis function `j`)". This is precisely what allows us to estimate the coefficients `coeffs_h_v[j]` for each basis function `j` as part of the HRF shape `h_v`.