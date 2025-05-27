# NEWS for hrfals

## 0.0.0.9000

* Added spatially regularised CF-ALS solver with graph Laplacian smoothing.
* New helper `build_voxel_laplacian()` constructs voxel neighbourhood Laplacians.
* Added `lambda_s`, `laplacian_obj` and solver parameters to `estimate_hrf_cfals` and `hrfals`.
* Convenience wrapper `estimate_hrf_spatial_cfals()` simplifies spatial smoothing usage.
* Benchmark function `benchmark_spatial_cfals()` explores performance for different penalties.

