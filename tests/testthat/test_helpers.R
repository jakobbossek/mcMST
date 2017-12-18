context("helper functions")

test_that("sampling of weights for weighted-sum approach", {
  weights = replicate(10L, sampleWeights(3L))
  expect_true(all(colSums(weights) == 1))
})

test_that("scalrization of weights matrices", {
  n = 10L
  x = matrix(10, nrow = n, ncol = n)
  y = matrix(10, nrow = n, ncol = n)
  z = matrix(100, nrow = n, ncol = n)

  expect_error(scalrizeWeights(list(x, y, z), c(0.3, 0.4)))
  expect_matrix(scalarizeWeights(list(x, y, z), sampleWeights(3L)), mode = "numeric", nrows = n, ncols = n)
  s = scalarizeWeights(list(x, y, z), c(0.9, 0.1, 0))
  expect_matrix(s, mode = "numeric", nrows = n, ncols = n)
  expect_true(all(s == 10))
})
