context("Solution enumerators")

test_that("Enumeration of Pareto-optimal solutions works", {
  n = 6L
  g = genRandomMCGP(n)

  res = getExactFront(g, obj.fun = objfunMCMST, enumerator.fun = enumerateMST, n.objectives = 2L)
  expect_matrix(res$pareto.front, nrows = 2L)
  expect_matrix(res$pareto.set, ncols = n - 2L)

  res = getExactFront(g, obj.fun = objfunMCTSP, enumerator.fun = enumerateTSP, n.objectives = 2L)
  expect_matrix(res$pareto.front, nrows = 2L)
  expect_matrix(res$pareto.set, ncols = n)
})
