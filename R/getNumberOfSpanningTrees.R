#' @title Compute number of spanning trees of a graph
#'
#' @description Makes use of Kirchhoff's matrix tree theorem to compute the
#' number of spanning trees of a given graph in polynomial time.
#'
#' @template arg_grapherator
#' @return [\code{integer(1)}]
#' @examples
#' # generate complete graph
#' g = genRandomMCGP(10)
#'
#' # this is equal to 10^8 (Cayley's theorem)
#' getNumberOfSpanningTrees(g)
#' @export
getNumberOfSpanningTrees = function(graph) {
  assertClass(graph, "grapherator")

  degrees = grapherator::getNodeDegrees(graph)
  # valence-matrix
  val.mat = diag(degrees)
  adj.mat = grapherator::getAdjacencyMatrix(graph)

  # laplace-matrix
  lap.mat = val.mat - adj.mat

  # drop arbitrary row and column
  lap.mat = lap.mat[-1L, -1L]

  n.strees = det(lap.mat)

  # due to numeric inaccuracy we get a double value
  n.strees = as.integer(ceiling(n.strees))
  return(n.strees)
}
