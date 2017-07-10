context("graph generators")

test_that("graph generation: simple 2o graph", {
  # here we generate a complex biobjective graph problem
  # with both euclidean and random weights

  g = mcGP(lower = 0, upper = 100)
  g = addCoordinates(g, n = 50L, generator = coordUniform)
  g = addWeights(g, method = "euclidean", symmetric = TRUE)
  g = addWeights(g, method = "random", weight.fun = runif, symmetric = TRUE)

  expect_class(g, "mcGP")
  expect_true(g$n.nodes == 50L)
  expect_true(g$n.clusters == 0L)
  expect_true(g$n.weights == 2L)
  expect_set_equal(g$weight.types, c("distance", "random"))
  expect_true(isSymmetricMatrix(g$weights[[1L]]))
  expect_true(isSymmetricMatrix(g$weights[[2L]]))
  expect_output(print(g), regexp = "MULTI")

  pls = plotGraph(g)
  expect_list(pls, types = "ggplot", len = 2L, any.missing = FALSE, all.missing = FALSE)
})

test_that("graph generation: complex clustered graph", {
  g = mcGP(lower = 0, upper = 100)
  g = addCenters(g, n.centers = 3L, generator = coordLHS)
  g = addCoordinates(g, n = c(5L, 10L, 15), by.centers = TRUE, generator = coordUniform, lower = c(0, 0), upper = c(1, 1))
  g = addCoordinates(g, n = 22, by.centers = TRUE, generator = coordUniform, lower = c(0, 0), upper = c(1, 1))
  g = addCoordinates(g, n = 100L, generator = coordGrid)
  g = addWeights(g, method = "random", weight.fun = rnorm, mean = 5, sd = 1.3)
  g = addWeights(g, method = "minkowski", p = 2.5, symmetric = FALSE)

  # check plotting of cluster centers
  pls = plotGraph(g, show.cluster.centers = TRUE)
  expect_list(pls, types = "ggplot", len = 2L, any.missing = FALSE, all.missing = FALSE)

  g = addWeights(g, method = "random", weight.fun = function(n) {
    sample(c(1, -10), n, replace = TRUE) * rexp(n, rate = 0.1) * 1:n
  })

  expect_class(g, "mcGP")
  expect_true(g$n.nodes == 152L)
  expect_true(g$n.clusters == 3L)
  expect_true(g$n.weights == 3L)
  expect_list(g$weights, types = "matrix", any.missing = FALSE, all.missing = FALSE, len = g$n.weights)
  expect_true(isSymmetricMatrix(g$weights[[1L]]))
  expect_true(isSymmetricMatrix(g$weights[[2L]]))
  expect_true(isSymmetricMatrix(g$weights[[3L]])) # distance based are always symmetric

  expect_error(plotGraph(g), regexpr = "not supported")
})

test_that("graph generation: manual passing of coordinates weights works", {
  g = mcGP(lower = 0, upper = 10)
  center.coordinates = matrix(c(1, 2, 2, 5, 8, 3), byrow = TRUE, ncol = 2L)
  g = addCenters(g, center.coordinates = center.coordinates)
  expect_equal(center.coordinates, g$center.coordinates)
  expect_true(g$n.clusters == nrow(center.coordinates))
  weights = diag(10)
  g = addWeights(g, weights = weights)
  g = addWeights(g, method = "random", weight.fun = rnorm, mean = 5, sd = 1.3)
  weights[1, 4] = 4
  g = addWeights(g, weights = weights)

  expect_class(g, "mcGP")
  expect_true(g$n.nodes == 10L)
  expect_true(g$n.clusters == 3L)
  expect_true(g$n.weights == 3L)
  expect_list(g$weights, types = "matrix", any.missing = FALSE, all.missing = FALSE, len = g$n.weights)
  expect_true(isSymmetricMatrix(g$weights[[1L]]))
  expect_true(isSymmetricMatrix(g$weights[[2L]]))
  expect_false(isSymmetricMatrix(g$weights[[3L]]))
})


test_that("graph generation: check correct error messages", {
  expect_error(mcGP(lower = 10, upper = 5))

  g = mcGP(lower = 0, upper = 100)
  expect_error(addWeights(g, method = "euclidean"), regexp = "number of nodes")
})
