context("Algorithms for the mcMST problem")

test_that("BG EMOA works well", {
  # bi-objective instance
  instance = genRandomMCGP(10L)

  res = mcMSTEmoaBG(instance, mu = 10L, max.iter = 50L)
  expect_class(res, "ecr_result")

  res = mcMSTEmoaBG(instance, mu = 10L, max.iter = 50L,
    mut = setup(mutEdgeExchange, instance = instance))
  expect_class(res, "ecr_result")

  # now increase number of weights
  instance = grapherator::addWeights(instance, generator = addWeightsRandom, method = rnorm, mean = 100, sd = 3)

  res = mcMSTEmoaBG(instance, mu = 10L, max.iter = 50L)
  expect_class(res, "ecr_result")

  res = mcMSTEmoaBG(instance, mu = 10L, max.iter = 50L,
    mut = setup(mutEdgeExchange, instance = instance))
  expect_class(res, "ecr_result")
})

test_that("Zhou EMOA works well", {
  instance = genRandomMCGP(10L)

  res = mcMSTEmoaZhou(instance, mu = 10L, ref.point = c(1e5, 1e5), max.iter = 50L)
  expect_class(res, "ecr_result")
})

test_that("Scalarization works well", {
  instance = genRandomMCGP(10L)

  res = mcMSTPrim(instance, n.lambdas = 10L)
  expect_list(res)
  expect_true(all(ecr::nondominated(res$pareto.front)))
})
