context("fastLSS benchmark")

# Simple test ensuring the benchmark helper works and
# fast implementation outperforms naive reference.

test_that("fastLSS is faster than naive implementation", {
  res <- benchmark_fastlss(n_seq = 20,
                           T_seq = 5,
                           v_seq = 3)
  expect_true(all(res$speedup > 1))
  expect_true(all(c("n", "T", "v", "fast", "naive", "memory", "speedup") %in%
                    names(res)))
})
