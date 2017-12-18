#' @title Metrics for spanning tree comparisson.
#'
#' @description Functions which expect two (spanning) trees and return a measure
#' of similiarity between those. Function \code{getNumberOfCommonEdges} returns
#' the (normalized) number of shared edges and function \code{getSizeOfLargestCommonSubtree}
#' returns the (normalized) size of the largest connected subtree which is located in
#' both trees.
#'
#' @param x [\code{matrix(2, n)}]\cr
#'   First spanning tree represented as a list of edges.
#' @param y [\code{matrix(2, n)}]\cr
#'   Second spanning tree represented as a list of edges.
#' @param n [\code{integer(1)} | \code{NULL}]\cr
#'   Number of nodes of the graph.
#'   Defaults to \code{length(x)}.
#' @param normalize [\code{logical(1)}]\cr
#'   Should measure be normalized to \eqn{[0, 1]} by devision
#'   through the number of edges?
#'   Default is \code{TRUE}.
#' @return [\code{numeric(1)}] Measure
#' @rdname similarity_metrics
#' @name similarity_metrics
#' @examples
#' # Here we generate two random spanning trees of a complete
#' # graph with 10 nodes
#' set.seed(1)
#' st1 = prueferToEdgeList(sample(1:10, size = 8, replace = TRUE))
#' st2 = prueferToEdgeList(sample(1:10, size = 8, replace = TRUE))
#' # Now check the number of common edges
#' NCE = getNumberOfCommonEdges(st1, st2)
#' # And the size of the largest common subtree
#' SLS = getSizeOfLargestCommonSubtree(st1, st2)
#' @export
getNumberOfCommonEdges = function(x, y, n = NULL, normalize = TRUE) {
  assertFlag(normalize)
  if (is.null(n))
    n = ncol(x) + 1L
  x.cv = edgeListToCharVec(x, n = n)
  y.cv = edgeListToCharVec(y, n = n)
  xy.cv = x.cv & y.cv
  if (normalize)
    return(sum(xy.cv / (n - 1L)))
  sum(xy.cv)
}

#' @rdname similarity_metrics
#' @export
getSizeOfLargestCommonSubtree = function(x, y, n = NULL, normalize = TRUE) {
  assertFlag(normalize)
  if (is.null(n))
    n = ncol(x) + 1L
  x.cv = edgeListToCharVec(x, n = n)
  y.cv = edgeListToCharVec(y, n = n)
  xy.cv = x.cv & y.cv
  xy = charVecToEdgelist(xy.cv)

  # now search for largest common, connected subtree
  nedges = ncol(xy)

  # here we get only the number and throw away the tree
  res = sapply(1:nedges, function(i) length(getReachableNodes(xy, xy[, i])))
  if (normalize)
    return((max(res) - 1L) / (n - 1))
  return(max(res) - 1L)
}

# Helper
#
# Given a edgelist and a node id, the function returns all nodes which are
# reachable starting at the given node.
getReachableNodes = function(edgelist, node) {
  nodes = node
  nodes2 = nodes
  while (length(nodes2) > 0L) {
    the.node = nodes2[1L]
    # check which edges are incident to current node
    is.adjacent = apply(edgelist, 2L, function(x) any(x == the.node))
    # get all unique nodes
    adjacent.nodes = unique(as.integer(edgelist[, is.adjacent]))
    nodes2 = setdiff(c(nodes2, adjacent.nodes), the.node)
    nodes = unique(c(nodes, adjacent.nodes))
    edgelist = edgelist[, !is.adjacent, drop = FALSE]
  }
  return(nodes)
}

