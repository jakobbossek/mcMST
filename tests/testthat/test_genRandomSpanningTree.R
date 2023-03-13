context("genRandomSpanningTrees")

test_that("random spanning tree generation returns reasonable results", {
  n = 10L
  tree = genRandomSpanningTree(n, type = "pruefer")
  expect_integer(tree, lower = 1, upper = n, len = n - 2L, any.missing = FALSE, all.missing = FALSE)

  tree = genRandomSpanningTree(n, type = "edgelist")
  expect_matrix(tree, mode = "integer", nrows = 2L, ncols = n - 1L, any.missing = FALSE, all.missing = FALSE)
  expect_true(all(tree >= 1 & tree <= n))

  tree = genRandomSpanningTree(n, type = "charvec")
  expect_integer(tree, lower = 0, upper = 1, len = n * n, any.missing = FALSE, all.missing = FALSE)
  expect_true(sum(tree) == (n - 1L))

  # Now multiple trees
  m = 3L
  trees.pruefer = genRandomSpanningTrees(m, n, type = "pruefer")
  expect_matrix(trees.pruefer, mode = "integer", nrows = m, ncols = n - 2L, any.missing = FALSE, all.missing = FALSE)
  trees.pruefer = genRandomSpanningTrees(m, n, type = "pruefer", simplify = FALSE)
  expect_list(trees.pruefer, types = "integer", len = m)

  trees.edgelist = genRandomSpanningTrees(m, n, type = "edgelist")
  expect_list(trees.edgelist, types = "matrix", len = m)
})
