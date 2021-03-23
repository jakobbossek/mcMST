grapheratorToGraph = function(g) {
  n = grapherator::getNumberOfNodes(g)
  w = grapherator::getNumberOfWeights(g)
  gpp = new(Graph, n, w, FALSE)
  for (i in 1:(n - 1L)) {
    for (j in (i + 1L):n) {
      if (g$adj.mat[i, j]) {
        ws = as.numeric(sapply(g$weights, function(w) w[i, j]))
        gpp$addEdge(i, j, ws, 0)
      }
    }
  }
  gpp$saveVectorOfEdges()
  gpp$setEdgeProbabilities(rep(1, gpp$getE()))
  return(gpp)
}
