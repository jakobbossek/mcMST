grapheratorToGraph = function(g) {
  n = grapherator::getNumberOfNodes(g)
  w = grapherator::getNumberOfWeights(g)
  weights = vector(mode = "numeric", length = w)
  gpp = new(Graph, n, w, FALSE)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j)
        next
      if (g$adj.mat[i, j]) {
        for (k in seq_len(w)) {
          weights[k] = as.numeric(g$weights[[k]][i, j])
        }
        gpp$addEdge(i, j, weights)
      }
    }
  }
  gpp$saveVectorOfEdges()
  gpp$setEdgeProbabilities(rep(1, gpp$getE()))
  return(gpp)
}
