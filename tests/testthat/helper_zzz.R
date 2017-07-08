set.seed(1)

isSymmetricMatrix = function(mat) {
  all(mat == t(mat))
}


objfunMCMST = function(pcode, instance) {
  getWeight(instance, prueferToEdgelist(pcode))
}

objfunMCTSP = function(perm, instance) {
  w1 = w2 = 0
  # add start node which is the end node as well
  perm = c(perm, 1L)
  # now add up weights
  for (i in 1:instance$n.nodes) {
    w1 = w1 + instance$weights[[1L]][perm[i], perm[i + 1L]]
    w2 = w2 + instance$weights[[2L]][perm[i], perm[i + 1L]]
  }
  return(c(w1, w2))
}
