context("crossover")

test_that("UnionCrossover produces valid solutions", {
  for (i in seq_len(10)) {
    g = genRandomMCGP(50L)
    gcpp = grapheratorToGraph(g)

    mst1 = gcpp$getRandomMST()
    mst2 = gcpp$getRandomMST()
    expect_true(mst1$isSpanningTree())
    expect_true(mst2$isSpanningTree())

    cx = ecr::setup(recUnionCrossover, instance = gcpp)

    mst = cx(list(mst1, mst2))
    expect_true(mst$isSpanningTree())

    max.degree = mst1$getMaximumDegree()

    cx = ecr::setup(recUnionCrossover, instance = gcpp, max.degree = max.degree)
    mst = cx(list(mst1, mst2))
    expect_true(mst$isSpanningTree())
    expect_true(mst$getMaximumDegree() <= max.degree)
  }
})
