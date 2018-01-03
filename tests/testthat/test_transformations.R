context("encoding transformations")

test_that("trafo of spanning tree encodings work well", {
  nl = 10; nu = 15
  for (i in seq_len(10L)) {
    n = sample(nl:nu, 1L)
    g = genRandomMCGP(n)
    pcode = genRandomPrueferCode(g)
    expect_length(pcode, n - 2L)

    # trafo to edgelist
    edgelist = prueferToEdgeList(pcode)
    expect_matrix(edgelist, nrows = 2L, ncols = n - 1L, mode = "numeric")
    expect_true(all(edgelist %in% 1:n))

    # convert to char vec and back
    cvec = edgeListToCharVec(edgelist, n = n)
    edgelist2 = charVecToEdgelist(cvec)

    # complicated way to check equality of edge lists
    # first sort each edge by node ids
    e1 = apply(edgelist, 2L, sort)
    e2 = apply(edgelist2, 2L, sort)
    # combine and search duplicates
    e12 = t(cbind(e1, e2))
    dups = duplicated(e12)
    # half of the edge must be duplicates
    expect_true(all(dups[n:(2 * n - 2)]))
    expect_true(all(!dups[1:(n-1)]))

    cvec1 = edgeListToCharVec(edgelist, n = n)
    cvec2 = prueferToCharVec(pcode)

    expect_true(all(cvec1 == cvec2))
    expect_true(all(cvec1 %in% c(0, 1)))
    expect_true(sum(cvec1) == n - 1L)
  }
})

test_that("trafo for roundtrip tours work well", {
  nl = 10; nu = 15
  for (i in seq_len(10L)) {
    n = sample(nl:nu, 1L)
    perm = genRandomPermutation(n)
    expect_length(perm, n)
    expect_true(all(sort(perm) == 1:n))

    edgelist = permutationToEdgelist(perm)
    expect_matrix(edgelist, nrows = 2L, ncols = n, mode = "numeric")
    expect_true(all(edgelist %in% 1:n))

    cvec1 = permutationToCharVec(perm, n = n)
    cvec2 = edgeListToCharVec(edgelist, n = n)

    expect_true(all(cvec1 == cvec2))
    expect_true(all(cvec1 %in% c(0, 1)))
  }
})

test_that("trafo between node list and edge list works", {
  for (i in 1:10L) {
    n = sample(2:10, 1L)
    nodelist = sample(1:10, size = n, replace = FALSE)
    edgelist = nodelistToEdgelist(nodelist)
    expect_matrix(edgelist, nrows = 2L, ncols = n - 1L, mode = "numeric")
    expect_true(all(as.integer(edgelist) %in% nodelist))
  }
})
