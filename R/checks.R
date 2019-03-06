# Function to check validity of subgraphs, spanning trees.

# Given a list of edge lists and a graph, this function checks whether all edges
# actually exist.
checkValidEdges = function(sols, g) {
  assertList(sols)
  assertClass(g, "grapherator")

  isValidSolution = function(sol, g) {
    #path = nodelistToEdgelist(path)
    edges.exist = apply(sol, 2L, function(edge) {
      g$adj.mat[edge[1L], edge[2L]]
    })
    all(edges.exist)
  }

  res = sapply(sols, isValidSolution, g)
  return(all(res))
}

# Check if each solution of a set is actually a spanning tree.
checkValidSpanningTrees = function(sols, g) {
  converter = new(RepresentationConverter)
  res = sapply(sols, function(sol) {
    tree = converter$edgeListToGraph(g, sol)
    tree$isSpanningTree()
  })
  return(invisible(all(res)))
}
