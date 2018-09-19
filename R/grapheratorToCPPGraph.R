grapheratorToCPPGraph = function(g) {
  n = grapherator::getNumberOfNodes(g)
  w = grapherator::getNumberOfWeights(g)
  gpp = new(Graph, n, w, FALSE)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j)
        next
      if (g$adj.mat[i, j]) {
        w1 = g$weights[[1]][i, j]
        w2 = g$weights[[2]][i, j]
        gpp$addEdge(i, j, w1, w2)
      }
    }
  }
  gpp$saveVectorOfEdges()
  return(gpp)
}
