context("build_voxel_laplacian")

library(Matrix)

# Simple 2-voxel mask 2x1x1
vol1 <- array(c(1, 1), dim = c(2, 1, 1))
res1 <- build_voxel_laplacian(vol1)

test_that("laplacian for two neighbouring voxels", {
  expect_equal(nrow(res1$L), 2)
  expect_equal(res1$degree, c(1, 1))
  expect_equal(as.matrix(res1$L), matrix(c(1, -1, -1, 1), 2, 2))
})

# 2x2x1 full mask
vol2 <- array(1, dim = c(2, 2, 1))
res2 <- build_voxel_laplacian(vol2)

test_that("degrees in 2x2x1 volume", {
  expect_equal(res2$degree, c(2, 2, 2, 2))
  expect_true(is(res2$L, "dgCMatrix"))
  expect_equal(sum(res2$L@x < 0), Matrix::nnzero(res2$L) / 2)
})

# Mask with isolated voxel
vol3 <- array(0, dim = c(2, 2, 1))
vol3[1,1,1] <- 1
vol3[2,2,1] <- 1
res3 <- build_voxel_laplacian(vol3)

test_that("isolated voxels have zero degree", {
  expect_equal(res3$degree, c(0, 0))
  expect_equal(Matrix::nnzero(res3$L), 0)
})
