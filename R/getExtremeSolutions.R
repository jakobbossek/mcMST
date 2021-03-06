#' @title Compute extreme spanning trees of bi-criteria graph problem.
#'
#' @description Internally \code{\link{mcMSTPrim}} is called with weights set accordingly.
#'
#' @template arg_grapherator
#' @return [\code{matrix(2, 2)}] The i-th column contains the objective vector
#' of the extreme i-th extreme solution
#' @export
getExtremeSolutions = function(graph) {
  assertClass(graph, "grapherator")
  if (grapherator::getNumberOfWeights(graph) != 2L)
    stopf("getExtremeSolutions: At the moment only bi-objective problems supported.")
  mcMSTPrim(graph, lambdas = c(0, 1))$pareto.front
}
