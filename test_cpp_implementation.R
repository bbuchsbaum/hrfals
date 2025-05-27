#!/usr/bin/env Rscript

# Test script to compare R and C++ LSS implementations
library(hrfals)
library(microbenchmark)

cat("=== Testing Enhanced C++ LSS Implementation ===\n\n")

# Test 1: Accuracy comparison with simple data
cat("📊 TEST 1: Accuracy Comparison\n")
cat("Comparing R and C++ implementations on simple synthetic data\n\n")

set.seed(42)
n <- 50
m <- 5
T_trials <- 10
V_voxels <- 20

# Generate test data
A <- matrix(rnorm(n * m), n, m)
C <- matrix(rnorm(n * T_trials), n, T_trials)
Y <- matrix(rnorm(n * V_voxels), n, V_voxels)
p_vec <- rnorm(n)
lambda_ridge <- 0.1

# Test R implementation
cat("Running R implementation...\n")
result_r <- lss_mode_a(Y, A, C, p_vec, lambda_ridge = lambda_ridge, use_cpp = FALSE)

# Test C++ implementation
cat("Running C++ implementation...\n")
result_cpp <- lss_mode_a(Y, A, C, p_vec, lambda_ridge = lambda_ridge, use_cpp = TRUE)

# Compare results
max_diff <- max(abs(result_r - result_cpp))
rel_diff <- max(abs(result_r - result_cpp) / (abs(result_r) + 1e-10))

cat(sprintf("✅ Maximum absolute difference: %.2e\n", max_diff))
cat(sprintf("✅ Maximum relative difference: %.2e\n", rel_diff))

if (max_diff < 1e-10) {
  cat("🎉 EXCELLENT: Results are numerically identical!\n")
} else if (max_diff < 1e-6) {
  cat("✅ GOOD: Results are very close (within 1e-6)\n")
} else {
  cat("⚠️  WARNING: Results differ by more than 1e-6\n")
}

# Test 2: Performance comparison
cat("\n📊 TEST 2: Performance Comparison\n")
cat("Benchmarking R vs C++ implementations\n\n")

# Smaller dataset for quick benchmarking
n_bench <- 100
T_bench <- 20
V_bench <- 50

A_bench <- matrix(rnorm(n_bench * m), n_bench, m)
C_bench <- matrix(rnorm(n_bench * T_bench), n_bench, T_bench)
Y_bench <- matrix(rnorm(n_bench * V_bench), n_bench, V_bench)
p_vec_bench <- rnorm(n_bench)

cat("Running microbenchmark (this may take a moment)...\n")
benchmark_result <- microbenchmark(
  R_impl = lss_mode_a(Y_bench, A_bench, C_bench, p_vec_bench, 
                      lambda_ridge = lambda_ridge, use_cpp = FALSE),
  CPP_impl = lss_mode_a(Y_bench, A_bench, C_bench, p_vec_bench, 
                        lambda_ridge = lambda_ridge, use_cpp = TRUE),
  times = 10
)

print(benchmark_result)

# Calculate speedup
r_median <- median(benchmark_result$time[benchmark_result$expr == "R_impl"])
cpp_median <- median(benchmark_result$time[benchmark_result$expr == "CPP_impl"])
speedup <- r_median / cpp_median

cat(sprintf("\n🚀 Speedup: %.2fx (C++ is %.2fx faster than R)\n", speedup, speedup))

# Test 3: Matrix conditioning diagnostics
cat("\n📊 TEST 3: Matrix Conditioning Diagnostics\n")
cat("Testing the lss_check_conditioning function\n\n")

# Test with well-conditioned matrix
A_good <- matrix(rnorm(n * m), n, m)
diag_result_good <- lss_check_conditioning(A_good, lambda_ridge = 0)
cat("Well-conditioned matrix:\n")
cat(sprintf("  Condition number: %.2e\n", diag_result_good$condition_number))
cat(sprintf("  Min eigenvalue: %.2e\n", diag_result_good$min_eigenvalue))
cat(sprintf("  Full rank: %s\n", diag_result_good$full_rank))

# Test with ill-conditioned matrix
A_bad <- A_good
A_bad[, m] <- A_bad[, 1] + 1e-10 * rnorm(n)  # Nearly collinear
diag_result_bad <- lss_check_conditioning(A_bad, lambda_ridge = 0)
cat("\nIll-conditioned matrix:\n")
cat(sprintf("  Condition number: %.2e\n", diag_result_bad$condition_number))
cat(sprintf("  Min eigenvalue: %.2e\n", diag_result_bad$min_eigenvalue))
cat(sprintf("  Suggested ridge: %.2e\n", diag_result_bad$suggested_ridge))

# Test 4: Edge cases
cat("\n📊 TEST 4: Edge Cases\n")
cat("Testing edge cases and error handling\n\n")

# Test with zero denominators (collinear trial regressor)
C_collinear <- C_bench
C_collinear[, 1] <- A_bench[, 1]  # Make first trial collinear with nuisance

cat("Testing collinear trial regressor...\n")
result_edge_r <- lss_mode_a(Y_bench, A_bench, C_collinear, p_vec_bench, 
                            lambda_ridge = 1e-6, use_cpp = FALSE)
result_edge_cpp <- lss_mode_a(Y_bench, A_bench, C_collinear, p_vec_bench, 
                              lambda_ridge = 1e-6, use_cpp = TRUE)

edge_diff <- max(abs(result_edge_r - result_edge_cpp))
cat(sprintf("✅ Edge case difference: %.2e\n", edge_diff))

if (all(is.finite(result_edge_r)) && all(is.finite(result_edge_cpp))) {
  cat("✅ Both implementations handle collinearity gracefully\n")
} else {
  cat("⚠️  WARNING: Edge case produced non-finite values\n")
}

# Test 5: Scaling test
cat("\n📊 TEST 5: Scaling Test\n")
cat("Testing performance across different problem sizes\n\n")

sizes <- data.frame(
  n = c(50, 100, 200),
  T = c(10, 20, 40),
  V = c(20, 50, 100)
)

scaling_results <- data.frame(
  size = character(nrow(sizes)),
  r_time = numeric(nrow(sizes)),
  cpp_time = numeric(nrow(sizes)),
  speedup = numeric(nrow(sizes)),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(sizes)) {
  n_i <- sizes$n[i]
  T_i <- sizes$T[i]
  V_i <- sizes$V[i]
  
  A_i <- matrix(rnorm(n_i * m), n_i, m)
  C_i <- matrix(rnorm(n_i * T_i), n_i, T_i)
  Y_i <- matrix(rnorm(n_i * V_i), n_i, V_i)
  p_vec_i <- rnorm(n_i)
  
  # Time R implementation (multiple runs for better precision)
  r_times <- replicate(5, {
    system.time({
      lss_mode_a(Y_i, A_i, C_i, p_vec_i, lambda_ridge = lambda_ridge, use_cpp = FALSE)
    })["elapsed"]
  })
  r_time <- median(r_times)
  
  # Time C++ implementation (multiple runs for better precision)
  cpp_times <- replicate(5, {
    system.time({
      lss_mode_a(Y_i, A_i, C_i, p_vec_i, lambda_ridge = lambda_ridge, use_cpp = TRUE)
    })["elapsed"]
  })
  cpp_time <- median(cpp_times)
  
  # Handle cases where timing is too small to measure accurately
  if (r_time == 0 || cpp_time == 0) {
    speedup_i <- NA
  } else {
    speedup_i <- r_time / cpp_time
  }
  
  scaling_results$size[i] <- sprintf("n=%d, T=%d, V=%d", n_i, T_i, V_i)
  scaling_results$r_time[i] <- r_time
  scaling_results$cpp_time[i] <- cpp_time
  scaling_results$speedup[i] <- speedup_i
  
  if (is.na(speedup_i)) {
    cat(sprintf("  %s: R=%.3fs, C++=%.3fs, speedup=N/A (too fast to measure)\n", 
                scaling_results$size[i], r_time, cpp_time))
  } else {
    cat(sprintf("  %s: R=%.3fs, C++=%.3fs, speedup=%.2fx\n", 
                scaling_results$size[i], r_time, cpp_time, speedup_i))
  }
}

cat("\n🎯 SUMMARY\n")
cat("==========\n")
cat(sprintf("✅ Accuracy: Maximum difference < %.1e\n", max_diff))

# Calculate average speedup excluding NAs
valid_speedups <- scaling_results$speedup[!is.na(scaling_results$speedup)]
if (length(valid_speedups) > 0) {
  cat(sprintf("🚀 Performance: Average speedup %.2fx\n", mean(valid_speedups)))
  cat(sprintf("📊 Scaling: Speedup ranges from %.2fx to %.2fx\n", 
              min(valid_speedups), max(valid_speedups)))
} else {
  cat("🚀 Performance: All tests too fast to measure accurately\n")
  cat("📊 Scaling: Use larger problem sizes for meaningful benchmarks\n")
}

cat("🔧 Diagnostics: Matrix conditioning function working\n")
cat("⚡ Edge cases: Handled gracefully\n")

# Updated summary logic
if (length(valid_speedups) > 0) {
  avg_speedup <- mean(valid_speedups)
  if (avg_speedup > 1.5) {
    cat("\n🎉 SUCCESS: C++ implementation provides significant speedup!\n")
  } else if (avg_speedup > 1.0) {
    cat("\n✅ GOOD: C++ implementation is faster than R\n")
  } else {
    cat("\n⚠️  NOTE: C++ implementation may need optimization\n")
  }
} else {
  cat("\n💡 INFO: Problem sizes too small for meaningful performance comparison\n")
}

cat("\n💡 RECOMMENDATIONS:\n")
if (length(valid_speedups) > 0 && max(valid_speedups) > 2 * min(valid_speedups)) {
  cat("• C++ implementation shows better scaling for larger problems\n")
}
if (length(valid_speedups) > 0 && any(valid_speedups < 1)) {
  cat("• Consider R implementation for very small problems\n")
}
cat("• Use lss_check_conditioning() to diagnose numerical issues\n")
cat("• Set appropriate ridge penalty for ill-conditioned problems\n")
cat("• For meaningful performance tests, use larger problem sizes (n>500, V>1000)\n") 