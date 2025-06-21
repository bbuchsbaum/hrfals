context("validate_hrf_engine_inputs")

# Use helper to create simple consistent data
simple <- simple_small_data()

# Successful validation returns correct dimensions

expect_dims <- validate_hrf_engine_inputs(simple$X_list, simple$Y,
                                          simple$Phi, simple$href)

test_that("validate_hrf_engine_inputs returns expected dims", {
  expect_equal(expect_dims$n, nrow(simple$Y))
  expect_equal(expect_dims$v, ncol(simple$Y))
  expect_equal(expect_dims$d, ncol(simple$X_list[[1]]))
  expect_equal(expect_dims$k, length(simple$X_list))
})

# Mismatched design matrix rows should error
X_bad <- simple$X_list
X_bad[[1]] <- rbind(X_bad[[1]], 0)


test_that("validate_hrf_engine_inputs detects row mismatch", {
  expect_error(
    validate_hrf_engine_inputs(X_bad, simple$Y, simple$Phi, simple$href),
    "Design matrices must have same rows as Y_proj"
  )
})

# Canonical reference must be normalized
href_bad <- simple$href * 2


test_that("validate_hrf_engine_inputs checks h_ref_shape_canonical normalization", {
  expect_error(
    validate_hrf_engine_inputs(simple$X_list, simple$Y, simple$Phi, href_bad),
    "must be normalised"
  )
})
