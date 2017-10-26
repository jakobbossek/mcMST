context("{write,read}GP")

test_that("writeGP and readGP work well", {
  # test instance
  g = mcGP(lower = c(0, 0), upper = c(100, 100))
  g = addCoordinates(g, n = 3, generator = coordUniform)
  g = addCoordinates(g, n = 9, by.centers = TRUE, generator = coordUniform, lower = c(0, 0), upper = c(10, 10))
  g = addEdges(g, method = "delauney")
  g = addWeights(g, method = "correlated", rho = -0.9)

  filename = "test.mcgp"
  writeGP(g, filename)
  g2 = readGP(filename)

  expect_equal(g$n.nodes, g2$n.nodes)
  expect_equal(g$n.clusters, g2$n.clusters)
  expect_equal(g$n.weights, g2$n.weights)
  expect_equal(g$membership, g2$membership)
  expect_equal(length(g$weights), length(g2$weights))
  expect_true(all(dim(g$coordinates) == dim(g2$coordinates)))
  unlink(filename)
})
