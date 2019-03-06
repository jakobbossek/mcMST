set.seed(1)

genRandomPrueferCode = function(g) {
  sample(1:g$n.nodes, size = g$n.nodes - 2L, replace = TRUE)
}

genRandomPermutation = function(n) {
  sample(1:n)
}

isSymmetricMatrix = function(mat) {
  all(mat == t(mat))
}

objfunMCMST = function(pcode, instance) {
  converter = new(RepresentationConverter)
  tree = converter$prueferCodeToGraph(instance, pcode)
  tree$getSumOfEdgeWeights()
}
