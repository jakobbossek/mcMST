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
  assertList(sols)
  assertClass(g, "grapherator")

  n = grapherator::getNumberOfNodes(g)

  # check existence of used edges
  if (!checkValidEdges(sols, g))
    stopf("There are solutions with non-existent edges!")

  for (sol in sols) {
    if (ncol(sol) != (n - 1))
      stopf("Each spanning tree of the graph must have %i edges, but solution has %i.", n - 1L, ncol(sol))

    if (!checkValidEdges(list(sol), g))
      stopf("There are solutions with non-existent edges!")

    # which nodes are reachable from node 1 (should be all
    # in a spanning tree)
    reachable.node.ids = getReachableNodes(sol, 1L)
    if (!setequal(1:n, reachable.node.ids)) {
      node.diff = setdiff(1:n, reachable.node.ids)
      stopf("Not all nodes in tree! The following nodes are not covered: %s", length(node.diff), collapse(node.diff, ", "))
    }
  }
  return(invisible(TRUE))
}
