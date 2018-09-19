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
getNumberOfCommonEdges = function(x, y, normalize = TRUE) {
  assertFlag(normalize)
  n.common = x$getNumberOfCommonEdges(y)

  if (!normalize)
    return(n.common)
  return(n.common / x$getE())
}

#' @rdname similarity_metrics
#' @export
getNumberOfCommonComponents = function(x, y, normalize = TRUE) {
  n.common = x$getNumberOfCommonComponents(y)
  return(n.common)
}

#' @rdname similarity_metrics
#' @export
getSizeOfLargestCommonSubtree = function(x, y, n = NULL, normalize = TRUE) {
  assertFlag(normalize)
  if (is.null(n))
    n = ncol(x) + 1L

  x = igraph::graph_from_edgelist(t(x), directed = FALSE)
  y = igraph::graph_from_edgelist(t(y), directed = FALSE)

  z = igraph::intersection(x, y)
  common.subtrees = igraph::components(z, mode = "weak")
  # subtract 1 since we are interested in number of edges
  sizes = common.subtrees$csize - 1L
  max.size = max(sizes)

  if (normalize)
    return(max.size / (n - 1))
  return(max.size)
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

#' @title Get common subtrees of two trees.
#'
#' @description Given two spanning trees, the function returns the subtrees
#' of the intersection of these.
#'
#' @param x [\code{matrix}]\cr
#'   Edge list of first tree.
#' @param y [\code{matrix}]\cr
#'   Edge list of second tree.
#' @param n [\code{integer(1)} | \code{NULL}]\cr
#'   Number of nodes.
#'   Default to \code{ncol(x) + 1}.
#' @return [\code{list}] List of matrizes. Each matrix contains the edges of one
#'   connected subtree.
#' @examples
#' # assume we have a graph with n = 10 nodes
#' n.nodes = 10
#' # we define two trees (matrices with colwise edges)
#' stree1 = matrix(c(1, 2, 1, 3, 2, 4, 5, 6, 6, 7), byrow = FALSE, nrow = 2)
#' stree2 = matrix(c(1, 3, 1, 2, 2, 4, 5, 8, 6, 7), byrow = FALSE, nrow = 2)
#' # ... and compute all common subtrees
#' subtrees = getCommonSubtrees(stree1, stree2, n = 10)
#' @export
getCommonSubtrees = function(x, y, n = NULL) {
  if (is.null(n))
    n = ncol(x) + 1L
  x.cv = edgeListToCharVec(x, n = n)
  y.cv = edgeListToCharVec(y, n = n)
  xy.cv = x.cv & y.cv
  xy = charVecToEdgelist(xy.cv)

  nedges = ncol(xy)

  # init component membership
  comps = rep(NA, nedges)
  comp = 1L
  for (i in 1:nedges) {
    # already in component
    if (!is.na(comps[i]))
      next
    nodes = xy[, i]
    reachable = getReachableNodes(xy, nodes)
    # which edges are used by the components?
    used.edges = apply(xy, 2L, function(edge) any(edge %in% reachable))
    comps[used.edges] = comp
    comp = comp + 1L
  }

  lapply(1:(comp - 1L), function(i) {
    xy[, which(comps == i), drop = FALSE]
  })
}
