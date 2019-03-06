context("Algorithms for the mcMST problem")

test_that("BG EMOA works well", {
  # bi-objective g
  g = genRandomMCGP(10L)
  gcpp = grapheratorToGraph(g)

  res = mcMSTEmoaBG(g, mu = 10L, max.iter = 50L)
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, gcpp)

  res = mcMSTEmoaBG(g, mu = 10L, max.iter = 50L,
    mut = mutKEdgeExchange)
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, gcpp)

  # now increase number of weights
  g = grapherator::addWeights(g, generator = addWeightsRandom, method = rnorm, mean = 100, sd = 3)
  gcpp = grapheratorToGraph(g)

  res = mcMSTEmoaBG(gcpp, mu = 10L, max.iter = 50L, mut = ecr::setup(mutSubgraphMST, scalarize = TRUE))
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, gcpp)

  res = mcMSTEmoaBG(gcpp, mu = 10L, max.iter = 50L,
    mut = mutKEdgeExchange)
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, gcpp)

  res = mcMSTEmoaBG(gcpp, mu = 10L, max.iter = 50L,
    mut = ecr::setup(mutKEdgeExchange, k = 2L, instance = gcpp))
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, gcpp)

  # check three objectives
  g = grapherator::addWeights(g, generator = addWeightsRandom, method = runif, min = 10, max = 100)
  gcpp = grapheratorToGraph(g)

  res = mcMSTEmoaBG(gcpp, mu = 10L, max.iter = 50L)
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, gcpp)
})

test_that("BG EMOA works on non-complete graphs", {
  g = grapherator::graph(0, 100)
  g = grapherator::addNodes(g, n = 5L, generator = grapherator::addNodesUniform)
  g = grapherator::addNodes(g, n = 4L, by.centers = TRUE, generator = grapherator::addNodesUniform, lower = c(0, 0), upper = c(10, 10))
  g = grapherator::addEdges(g, type = "intracluster", generator = grapherator::addEdgesDelauney)
  g = grapherator::addEdges(g, type = "intercluster", generator = grapherator::addEdgesSpanningTree, runs = 2L, k = 3L)
  g = grapherator::addWeights(g, generator = grapherator::addWeightsCorrelated, rho = -0.6)
  gcpp = grapheratorToGraph(g)

  res = mcMSTEmoaBG(gcpp, mu = 10L, max.iter = 50L, mut = ecr::setup(mutSubgraphMST, instance = gcpp, scalarize = TRUE))
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, gcpp)

  res = mcMSTEmoaBG(gcpp, mu = 10L, max.iter = 50L, mut = ecr::setup(mutKEdgeExchange, instance = gcpp))
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, gcpp)
})

test_that("Zhou EMOA works well", {
  g = genRandomMCGP(10L)
  gcpp = grapheratorToGraph(g)

  res = mcMSTEmoaZhou(gcpp, mu = 10L, ref.point = c(1e5, 1e5), max.iter = 50L)
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, gcpp)
})

test_that("Scalarization works well", {
  g = genRandomMCGP(25L)
  gcpp = grapheratorToGraph(g)

  res = mcMSTWeightedSum(gcpp, n.lambdas = 10L)
  expect_list(res)
  expect_true(all(ecr::nondominated(t(as.matrix(res$pareto.front)))))
  checkValidSpanningTrees(res$pareto.set, gcpp)
})
