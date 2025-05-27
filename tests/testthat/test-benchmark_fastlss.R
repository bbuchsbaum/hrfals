context("fastLSS benchmark")

# Test ensuring the benchmark helper function works correctly
# and produces valid timing and memory usage results.

test_that("fastLSS benchmark function works correctly", {
  res <- benchmark_fastlss(n_seq = 50,
                           T_seq = 20,
                           v_seq = 10)
  
  # Check that the benchmark function returns the expected structure
  expect_true(all(c("n", "T", "v", "fast", "naive", "memory", "speedup") %in%
                    names(res)))
  
  # Check that we have the expected number of rows
  expect_equal(nrow(res), 1)
  
  # Check that timing values are non-negative
  expect_true(all(res$fast >= 0))
  expect_true(all(res$naive >= 0))
  
  # Check that memory usage is positive
  expect_true(all(res$memory > 0))
  
  # Check that speedup is computed correctly (allowing for numerical precision)
  expected_speedup <- res$naive / res$fast
  expect_equal(res$speedup, expected_speedup, tolerance = 1e-10)
})

test_that("fastLSS benchmark works with multiple problem sizes", {
  # Test multiple problem sizes to ensure benchmark function handles variety
  res <- benchmark_fastlss(n_seq = c(100, 200),
                           T_seq = c(50, 100),
                           v_seq = c(20))
  
  # Check that we have results for all combinations
  expect_equal(nrow(res), 4)
  
  # Check that all required columns are present
  expect_true(all(c("n", "T", "v", "fast", "naive", "memory", "speedup") %in%
                    names(res)))
  
  # Check that timings are non-negative
  expect_true(all(res$fast >= 0))
  expect_true(all(res$naive >= 0))
  
  # Check that memory usage is positive
  expect_true(all(res$memory > 0))
  
  # Check that speedup calculation is correct
  expected_speedups <- res$naive / res$fast
  expect_equal(res$speedup, expected_speedups, tolerance = 1e-10)
  
  # Check that speedup values are finite where both timings are positive
  # (handle cases where timing might be 0.000 due to measurement precision)
  valid_timings <- res$fast > 0 & res$naive > 0
  if (any(valid_timings)) {
    expect_true(all(is.finite(res$speedup[valid_timings])))
  }
})
