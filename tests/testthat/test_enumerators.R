context("Solution enumerators")

test_that("Enumeration of Pareto-optimal solutions works", {
  n = 6L
  g = genRandomMCGP(n)
  gcpp = grapheratorToGraph(g)

  # check on complete graph
  res = getExactFront(gcpp, obj.fun = objfunMCMST, enumerator.fun = enumerateMST, n.objectives = 2L)
  expect_data_frame(res$pareto.front, ncols = 2L)
  expect_matrix(res$pareto.set, ncols = n - 2L)

  # check on degenerated graph with only one spanning tree
  #FIXME: reenable
  # g = grapherator::graph(0, 100)
  # g = grapherator::addNodes(g, n = n, generator = grapherator::addNodesUniform)
  # g = grapherator::addEdges(g, generator = grapherator::addEdgesSpanningTree)
  # g = grapherator::addWeights(g, generator = grapherator::addWeightsRandom, method = runif, min = 10, max = 100)
  # g = grapherator::addWeights(g, generator = grapherator::addWeightsRandom, method = runif, min = 10, max = 100)
  # gcpp = grapheratorToGraph(g)

  # res = getExactFront(gcpp, obj.fun = objfunMCMST, enumerator.fun = enumerateMST, n.objectives = 2L, simplify = FALSE)
  # expect_matrix(res$pareto.front, nrows = 2L, ncols = 1L)
  # expect_matrix(res$pareto.set, nrows = 2L)
})
