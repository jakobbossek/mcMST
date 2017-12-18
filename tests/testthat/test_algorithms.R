context("Algorithms for the mcMST problem")

test_that("BG EMOA works well", {
  # bi-objective g
  g = genRandomMCGP(10L)

  res = mcMSTEmoaBG(g, mu = 10L, max.iter = 50L)
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, g)

  res = mcMSTEmoaBG(g, mu = 10L, max.iter = 50L,
    mut = setup(mutEdgeExchange, instance = g))
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, g)

  # now increase number of weights
  g = grapherator::addWeights(g, generator = addWeightsRandom, method = rnorm, mean = 100, sd = 3)

  res = mcMSTEmoaBG(g, mu = 10L, max.iter = 50L, scalarize = TRUE)
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, g)

  res = mcMSTEmoaBG(g, mu = 10L, max.iter = 50L,
    mut = setup(mutEdgeExchange, instance = g))
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, g)

  # check three objectives
  g = grapherator::addWeights(g, generator = addWeightsRandom, method = runif, min = 10, max = 100)
  res = mcMSTEmoaBG(g, mu = 10L, max.iter = 50L)
  expect_class(res, "ecr_result")
  checkValidSpanningTrees(res$pareto.set, g)
})

test_that("Zhou EMOA works well", {
  g = genRandomMCGP(10L)

  res = mcMSTEmoaZhou(g, mu = 10L, ref.point = c(1e5, 1e5), max.iter = 50L)
  expect_class(res, "ecr_result")
  res$pareto.set = lapply(res$pareto.set, prueferToEdgeList)
  checkValidSpanningTrees(res$pareto.set, g)
})

test_that("Scalarization works well", {
  g = genRandomMCGP(10L)

  res = mcMSTPrim(g, n.lambdas = 10L)
  expect_list(res)
  expect_true(all(ecr::nondominated(res$pareto.front)))
})
