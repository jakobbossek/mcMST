context("EMOAs for mcMST problem")

test_that("BG EMOA works well", {
  instance = genRandomMCGP(10L)

  res = emoaMST_BG(instance, mu = 10L, ref.point = c(1e5, 1e5), max.iter = 50L)
  expect_class(res, "ecr_result")
})

test_that("Zhou EMOA works well", {
  instance = genRandomMCGP(10L)

  res = emoaMST_BG(instance, mu = 10L, ref.point = c(1e5, 1e5), max.iter = 50L)
  expect_class(res, "ecr_result")
})
