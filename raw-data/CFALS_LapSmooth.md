Okay, this is very helpful "meta-review" guidance. I'll incorporate these final suggestions, especially the API strategy (Option A+B) and the refined implementation order, into the final proposal document.

---

**Final Proposal: Voxel-wise Spatial Smoothness Regularization for CF-ALS**

**1. Introduction and Motivation**

(Retain from previous version)

The Confound-Free Alternating Least Squares (CF-ALS) method, as implemented in the `hrfals` package, provides a robust approach for data-driven estimation of Hemodynamic Response Functions (HRFs) from fMRI data. It estimates HRF basis coefficients (`H`) and condition-specific amplitudes (`B`) by iteratively solving regularized least-squares problems.

Currently, CF-ALS estimates HRF coefficients for each voxel independently. While effective, this can lead to noisy and spatially implausible HRF estimates, especially in data with moderate-to-high noise levels or when using flexible HRF bases (e.g., FIR or B-splines with many knots). Adjacent voxels, which often share similar physiological responses, may exhibit divergent HRF shapes.

This proposal outlines the addition of a spatial smoothness penalty to the CF-ALS objective function. This penalty will encourage the HRF coefficient vectors of neighboring voxels to be similar, leveraging the inherent spatial coherence of fMRI data. The specific form proposed is a graph Laplacian regularizer.

The benefits of this addition are:
*   **Reduced Variance:** Spatial regularization will reduce the variance of HRF estimates, leading to more stable and reliable results.
*   **Improved Anatomical Plausibility:** Smoother HRF maps are generally more consistent with the underlying neurophysiology.
*   **Preservation of Boundaries:** The Laplacian penalty, particularly the unnormalized version, is known to preserve sharp functional boundaries when the regularization strength (`λs`) is appropriately tuned, more so than simple Gaussian smoothing.
*   **Backward Compatibility:** The proposed changes will leave the current CF-ALS behavior untouched when the spatial regularization strength `λs` is set to zero.
*   **Single New Hyperparameter:** Introduces only one new user-tunable hyperparameter, `λs`.

**2. Current CF-ALS Objective and h-update (Relevant Parts)**

The existing CF-ALS algorithm (as seen in `R/cf_als_engine.R`) minimizes an objective function. For the h-update, the parameters `lambda_h` (for `R_eff` or identity penalty) and `lambda_joint` (for an additional identity penalty) are used. `lambda_joint` remains optional and defaults to 0; we recommend `lambda_joint ≥ 0.1` when using a high-order FIR basis (many basis functions) to help keep the per-voxel `LHS_v` matrix well-conditioned, especially if `lambda_h` is small.

The h-update step in `cf_als_engine.R` for each voxel `v` (column `H_v` of `H`) solves:
`LHS_v H_v = RHS_v`
where:
*   `LHS_v = lambda_h * h_penalty_matrix + lambda_joint * diag(d) + sum_l (B_{lv}^2 * X_l^T X_l)` (if `!fullXtX_flag`).
    *   `h_penalty_matrix` is `R_mat_eff` (if provided) or `diag(d)` (if `R_mat_eff` is NULL).
*   `RHS_v = sum_l (B_{lv} * X_l^T Y_v)`.

This update is performed independently for each voxel.

**3. Proposed Spatially Regularized Objective**

We propose to add a spatial penalty term to the global objective function:

`min_{H,B} Σ_v ||Y_v - Σ_c (X_c H_v) B_{cv}||_2^2 + λ_b||B||_F^2 + λ_h||H R_{eff}^{1/2}||_F^2 + λ_joint||H||_F^2 + λ_s tr(H L H^T)`

where:
*   `H` is the `d x v` matrix of HRF basis coefficients.
*   `L` is the `v x v` **unnormalized** graph Laplacian: `L = Deg - Adj`.
*   `λ_s` is the non-negative scalar for spatial regularization.

**4. Derivation of the New h-update**

The β-update remains unchanged.
The spatially regularized h-update couples all voxels. For `vec(H)` (stacking columns of `H`):

Let `A_H = BlockDiag(LHS_v) + λ_s (L ⊗ I_d)`. The system is `A_H vec(H) = vec(RHS_matrix)`.
Here:
*   `RHS_matrix` is `d x v`, with column `v` being `RHS_v`.
*   `BlockDiag(LHS_v)` has `d x d` blocks `LHS_v` (including `λ_h R_eff` and `λ_joint I_d`) on its diagonal.
*   `(L ⊗ I_d)` is the Kronecker product, applying spatial coupling.

**5. Solving the Coupled h-update System**

`A_H` is symmetric and positive definite for `λ_h > 0, λ_s ≥ 0` (assuming `R_eff` is at least PSD).

**Solver Strategies:**

1.  **Small `v*d` (e.g., `< 50,000`): Direct Sparse Solver**
    *   Construct `A_H = BlockDiag(list_LHS_v_blocks) + lambda_s * kronecker(L_mat, diag(d))` explicitly as a sparse `dgCMatrix`.
    *   Solve using `Matrix::solve(A_H, as.vector(RHS_matrix))`.

2.  **Large `v*d`: Conjugate Gradient (CG) Method using `Rlinsolve`**
    *   The `Rlinsolve` package provides efficient iterative solvers for sparse linear systems, including `lsolve.cg`.
    *   **Matrix Construction:** Construct `A_H` as above (sparse `dgCMatrix`).
    *   **Invocation:** Call `Rlinsolve::lsolve.cg(A = A_H, B = as.vector(RHS_matrix), preconditioner = P_sparse, tol = cg_tol, maxiter = cg_max_iter, xinit = as.vector(h_current), adjsym = TRUE, verbose = FALSE)`. The `adjsym = TRUE` is appropriate as `A_H` is symmetric.
    *   **Preconditioner `P_sparse`:**
        *   Construct a sparse block-diagonal preconditioner matrix `P`. Each `d x d` diagonal block `P_v` would be `solve(LHS_vx_block + lambda_s * laplacian_obj$degree[vx] * diag(d))`.
        *   `P_sparse = BlockDiag(list_of_P_v_blocks)`. This `P_sparse` matrix is then passed to `lsolve.cg`.
    *   **CG Tolerance & Iterations:** `cg_tol = 1e-4`. `cg_max_iter` (e.g., 100). If outer ALS sweep parameter change is `< 1e-3`, cap `cg_max_iter` at a lower value (e.g., 20).

**6. Implementation Details**

*   **New Helper Function: `build_voxel_laplacian(volume, connectivity = 6)`**
    *   Input: `volume` (a `neuroim2::NeuroVol` object, typically a `LogicalNeuroVol` mask).
    *   Output: `list(L = L_sparse_matrix, degree = degree_vector_numeric_v)`.
    *   Handles mapping between 3D and 1D indices. Robust to isolated/disconnected voxels in the mask (e.g., by ensuring `degree` entries are non-zero for included voxels, perhaps by adding a small epsilon to `L`'s diagonal if necessary for solver stability, or by filtering out zero-degree voxels from the analysis if they are not truly part of the intended spatial domain).

*   **New Internal Helper: `make_lhs_block_list(params...) -> list_of_dxd_matrices`**
    *   Extracts the logic for constructing `LHS_v` (the `d x d` matrix for each voxel before spatial coupling) into a reusable function. This list of matrices is computed once per h-update sweep in an ALS iteration.

*   **Modifications to `R/cf_als_engine.R` (`cf_als_engine`)**
    *   Add `Rlinsolve` to `Imports` in `DESCRIPTION`.
    *   Add parameters: `lambda_s = 0`, `laplacian_obj = NULL`, `cg_max_iter = 100`, `cg_tol = 1e-4`, `h_solver = "auto"`, `...`.
    *   Inside `cf_als_engine`, after β-update:
        ```R
        if (lambda_s == 0 || is.null(laplacian_obj) || is.null(laplacian_obj$L)) {
            # <Current per-voxel h-update loop, using make_lhs_block_list() internally for each LHS_v>
        } else {
            L_mat <- laplacian_obj$L
            degree_vec <- laplacian_obj$degree 

            list_LHS_v_blocks <- make_lhs_block_list(X_list_proj, b_current, lambda_h, lambda_joint, R_mat_eff, fullXtX_flag, ...) 
            RHS_matrix <- # ... calculate (d x v) RHS_matrix ...

            # Construct A_H (sparse dgCMatrix)
            # A_H_spatial <- construct_A_H_sparse(list_LHS_v_blocks, lambda_s, L_mat, d, v)
            # (Helper function to build A_H = BlockDiag(list_LHS_v_blocks) + lambda_s * kronecker(L_mat, diag(d)))

            current_solver <- h_solver
            if (current_solver == "auto") {
                # Heuristic: if A_H would have > ~30M non-zeros, consider CG, else direct.
                # nnz_A_H_est <- v*d^2 + Matrix::nnzero(L_mat)*d 
                # current_solver <- if (nnz_A_H_est < 30e6 && d * v < 150000) "direct" else "cg" 
                 current_solver <- if (d * v < 50000) "direct" else "cg" # Simpler heuristic for now
            }

            if (current_solver == "direct") {
                # vec_H_new = as.vector(Matrix::solve(A_H_spatial, as.vector(RHS_matrix)))
            } else { # current_solver == "cg"
                # Construct sparse preconditioner matrix P_sparse
                # list_P_v_blocks <- lapply(1:v, function(vx_idx) {
                #    LHS_diag_block <- list_LHS_v_blocks[[vx_idx]]
                #    P_block <- solve(LHS_diag_block + lambda_s * degree_vec[vx_idx] * diag(d))
                #    return(P_block)
                # })
                # P_sparse <- Matrix::bdiag(list_P_v_blocks) # Forms a sparse block-diagonal matrix
                
                # cg_sol <- Rlinsolve::lsolve.cg(A = A_H_spatial, 
                #                               B = as.vector(RHS_matrix), 
                #                               preconditioner = P_sparse, # May need to be NULL if P_sparse is too costly or ineffective
                #                               xinit = as.vector(h_current), 
                #                               reltol = cg_tol, # Rlinsolve uses 'reltol'
                #                               maxiter = cg_max_iter, 
                #                               adjsym = TRUE, verbose = FALSE) 
                # vec_H_new = as.vector(cg_sol$x)
                # if (cg_sol$iter == cg_max_iter && cg_sol$errors[length(cg_sol$errors)] > cg_tol) {
                #    warning("CG solver did not converge within max_iter for h-update.")
                # }
            }
            h_current = matrix(vec_H_new, d, v)
        }
        # Apply per-voxel normalization and sign alignment (normalize_and_align_hrf)
        ```

*   **API Strategy (Option A+B):**
    1.  Modify `estimate_hrf_cfals` (and its internal call to `cf_als_engine`) to include the new parameters (`lambda_s`, `laplacian_obj`, etc.). This is the primary, unified engine.
    2.  Add a new exported convenience wrapper:
        `estimate_hrf_spatial_cfals(..., lambda_s_default = 0.1, laplacian_obj, ...)`
        This wrapper will simply call `estimate_hrf_cfals` with the provided `laplacian_obj` and a sensible non-zero default for `lambda_s`, forwarding other arguments.

*   **Modifications to Wrapper Functions (`R/hrfals_control_defaults.R`, etc.)**:
    *   Update defaults, signatures, and pass `...`.

**7. Parameter `λs`**

*   Users will typically create `laplacian_obj` via `build_voxel_laplacian(neuroim_mask_vol)` and pass it.
*   Guidance on `λs` selection (cross-validation, visual inspection).
*   API documentation for `lambda_s` should note: "The effective physical smoothness depends on `lambda_s` relative to voxel size. For comparability across resolutions, consider scaling `lambda_s` inversely with squared voxel edge length (e.g., `lambda_s_phys = lambda_s / mean(voxel_edge_length)^2`)."

**8. Impact on Existing Code & S3 Methods**

(Retain from previous version, still valid)

**9. Unit Tests**

(Retain from previous version, with note on `neuroim2::LogicalNeuroVol` for `build_voxel_laplacian`)

**10. Computational Impact**

*   **Direct solver memory (worst-case):** `A_H` construction involves `v` dense `d x d` blocks and `d` copies of `L`. `nnz(A_H) ≈ v*d^2 + nnz(L)*d`.
*   **CG MVP cost:** Dominated by `v*d^2` (block-diagonal part using pre-stored `LHS_v` blocks) and `d*nnz(L)` (Laplacian part). Real wall-time often memory-bandwidth bound.

**11. Risks and Mitigations**

(Retain from previous version, with updated `Deg_{vv}` detail for preconditioner)

**12. Implementation Order & Ticketed Sprint Plan**

This plan aims for incremental integration and testing.

---

**Sprint Plan: Spatial CF-ALS Integration**

**Total Estimated Time: 4 Weeks**

**Week 0.5: Preparatory Refactoring (MR #0)**
*   **Ticket HALS-S01:** Refactor `LHS_v` Block Construction.
    *   **Task:** Create internal helper `make_lhs_block_list(X_list_proj, b_current, lambda_h, lambda_joint, R_mat_eff, fullXtX_flag, d, k, v, ...)` in `R/cf_als_engine.R` (or a utils file) that computes and returns a list of all `v` (d x d) `LHS_v` matrices.
    *   **Task:** Modify the existing per-voxel h-update loop in `cf_als_engine.R` to call `make_lhs_block_list()` once per h-sweep and then use the precomputed blocks.
    *   **Verification:** Ensure all existing tests for `cf_als_engine` pass with this refactoring (no change in output).
    *   **Deliverable:** MR with refactored `cf_als_engine.R`.

**Week 1: Laplacian Builder & Initial CG Integration (MR #1)**
*   **Ticket HALS-S02:** Implement `build_voxel_laplacian`.
    *   **Task:** Create `R/utils_spatial.R` (or similar). Implement `build_voxel_laplacian(volume, connectivity=6)` taking a `neuroim2::NeuroVol` (mask) and returning `list(L = dgCMatrix, degree = numeric_vector_v)`.
    *   **Task:** Add unit tests for `build_voxel_laplacian` with small synthetic `LogicalNeuroVol` masks (e.g., 2x1x1, 2x2x1, handling of isolated voxels).
*   **Ticket HALS-S03 (Revised):** Integrate `Rlinsolve::lsolve.cg` for h-update.
    *   **Task:** Add `Rlinsolve` to package Imports.
    *   **Task:** Implement the logic to construct the sparse `A_H` matrix.
    *   **Task:** Implement the logic to construct the sparse block-Jacobi preconditioner matrix `P_sparse`.
    *   **Task:** Add the branch in `cf_als_engine.R` to use `Rlinsolve::lsolve.cg` with `A_H` and `P_sparse` when `lambda_s > 0 && h_solver == "cg"`.
    *   **Task:** Basic numerical verification on a small (e.g., 20-200 voxel) toy problem, comparing against the direct solver path if also implemented or a known solution.
    *   **Deliverable:** MR with `build_voxel_laplacian` (from HALS-S02) and `Rlinsolve` integration for CG path (from HALS-S03 Revised).

**Week 2: Direct Solver, Wrappers, Alias (MR #2)**
*   **Ticket HALS-S04 (Combined into S05):** (No separate Rcpp CG implementation needed if `Rlinsolve` is used).
*   **Ticket HALS-S05 (Revised):** Implement Direct Solver Path & Finalize CG path.
    *   **Task:** Implement the "direct" solver path in `cf_als_engine.R` (construct sparse `A_H` and use `Matrix::solve()`).
    *   **Task:** Refine `h_solver="auto"` logic based on initial performance tests of `A_H` construction and solve times.
    *   **Task:** Ensure CG and Direct solver paths are robust.
*   **Ticket HALS-S06:** Update Wrappers and Add Alias.
    *   **Task:** Add `lambda_s`, `laplacian_obj`, `h_solver`, `cg_max_iter`, `cg_tol`, `...` to `hrfals_control_defaults`, `estimate_hrf_cfals`, and `hrfals` function signatures.
    *   **Task:** Implement the new convenience wrapper `estimate_hrf_spatial_cfals(...)` that calls `estimate_hrf_cfals` with a default `lambda_s > 0`.
    *   **Task:** Add basic unit tests for parameter pass-through and alias functionality.
    *   **Deliverable:** MR with Direct Solver Path, finalized CG path (from HALS-S05 Revised), and updated wrappers/alias (from HALS-S06).

**Week 3: Comprehensive Testing & Benchmarking (MR #3)**
*   **Ticket HALS-S07:** Extended Unit & Integration Tests.
    *   **Task:** Implement all unit tests outlined in Proposal Section 9 (numerical identity, small examples, opposite sign checks, solver comparisons).
    *   **Task:** Test with `lambda_s=0` to ensure exact match with non-spatial version.
*   **Ticket HALS-S08:** Performance Benchmarking.
    *   **Task:** Adapt `test-method_comparison_revised.R` to include `lambda_s` variations.
    *   **Task:** Profile direct vs. CG solver performance for different `v*d` sizes.
    *   **Task:** Benchmark on 1-2 real (or realistic large synthetic) datasets.
    *   **Deliverable:** MR with comprehensive tests and an internal report/plots summarizing benchmark results and effect of `lambda_s`.

**Week 4: Documentation & Finalization (MR #4)**
*   **Ticket HALS-S09:** Update Documentation.
    *   **Task:** Update Roxygen documentation for all modified/new functions and parameters (especially `lambda_s` scaling note, `laplacian_obj` usage).
    *   **Task:** Update `NEWS.md`.
*   **Ticket HALS-S10:** Create Vignette/Example.
    *   **Task:** Add a section to an existing vignette or create a small new one demonstrating "Spatially-Smoothed CF-ALS". Include an example of `build_voxel_laplacian(neuroim_mask_vol)` and passing `laplacian_obj` to `estimate_hrf_spatial_cfals`.
    *   **Task:** (Optional, if time permits) Add a small example for auto-tuning `lambda_s` on a subset of data or using a heuristic.
    *   **Deliverable:** MR with final documentation and vignette.