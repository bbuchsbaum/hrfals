context("auto_residualize")

set.seed(123)

n <- 30
Tt <- 5

# Case m small -> Woodbury
m_small <- 4
A_small <- matrix(rnorm(n * m_small), n, m_small)
C <- matrix(rnorm(n * Tt), n, Tt)

V_ref_small <- woodbury_residualize(C, A_small)
V_auto_small <- auto_residualize(C, A_small, woodbury_thresh = 50)

test_that("auto_residualize uses Woodbury below threshold", {
  expect_equal(V_auto_small, V_ref_small, tolerance = 1e-12)
})

# Case m large -> QR
m_large <- 60
A_large <- matrix(rnorm(n * m_large), n, m_large)
V_ref_large <- qr_residualize(C, A_large)
V_auto_large <- auto_residualize(C, A_large, woodbury_thresh = 50)

test_that("auto_residualize falls back to QR above threshold", {
  expect_equal(V_auto_large, V_ref_large, tolerance = 1e-12)
})
