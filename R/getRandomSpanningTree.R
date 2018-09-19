#' @title Generate random spanning tree.
#'
#' @description Given a \code{\link[grapherator]{grapherator}} object this function
#' returns a random spanning tree. The tree generation process is a simple heuristic:
#' A random weight from a \eqn{U(0, 1)}-distribution is assigned to each edge of the
#' graph. Next, a spanning tree is computed by \code{\link[vegan]{spantree}}.
#'
#' @note Most likely this heuristic does not produce each spanning tree with equal
#' probability.
#'
#' @template arg_grapherator
#' @return [\code{matrix}] Edge list of spanning tree edges.
#' @examples
#' g = genRandomMCGP(10L)
#' stree = getRandomSpanningTree(g)
#' @export
getRandomSpanningTree = function(graph) {
  graph$getMST()
}
