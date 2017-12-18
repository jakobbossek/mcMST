context("similarity measures")

test_that("similarity measures work", {
  n.nodes = 10L

  getRandomPCode = function(n) {
    sample(1:n, size = n - 2L, replace = TRUE)
  }
  # check on random spanning trees
  for (i in 1:n.nodes) {
    stree1 = prueferToEdgeList(getRandomPCode(n.nodes))
    stree2 = prueferToEdgeList(getRandomPCode(n.nodes))
    expect_number(getNumberOfCommonEdges(stree1, stree2), lower = 0, upper = 1)
  }

  # build trees by hand
  # Edges 1, 2, 3 and 5 are in both trees
  # First three edges form the largest connected component
  stree1 = matrix(c(1, 2, 1, 3, 2, 4, 5, 6, 6, 7), byrow = FALSE, nrow = 2L)
  stree2 = matrix(c(1, 3, 1, 2, 2, 4, 5, 8, 6, 7), byrow = FALSE, nrow = 2L)

  NCE = getNumberOfCommonEdges(stree1, stree2, n = n.nodes, normalize = FALSE)
  expect_number(NCE, lower = 0, upper = n.nodes)
  expect_true(NCE == 4)

  SLS = getSizeOfLargestCommonSubtree(stree1, stree2, n = n.nodes, normalize = FALSE)
  expect_number(SLS, lower = 0, upper = n.nodes)
  expect_true(SLS == 3)
})

test_that("similarity matrix calculation works", {
  n = 10L
  mu = 10L
  g = genRandomMCGP(n)

  # approximate front
  res = mcMSTEmoaBG(g, mu = mu, max.iter = 50L)
  set = res$pareto.set

  sim.mat = computeSimilarityMatrix(set, sim.fun = getNumberOfCommonEdges, normalize = FALSE)
  expect_matrix(sim.mat, mode = "numeric", nrows = mu, ncols = mu, any.missing = FALSE, all.missing = FALSE)
  # sim(i, i) = n - 1
  expect_true(all(diag(sim.mat) == (n - 1L)))
  # sim(i, j) = sim(j, i)
  expect_true(isSymmetric(sim.mat))
})
