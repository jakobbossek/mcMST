#' @title Multi-Objective Prim algorithm.
#'
#' @description Approximates the Pareto-optimal mcMST front of a multi-objective
#' graph problem by iteratively applying Prim's algorithm for the single-objective
#' MST problem to a scalarized version of the problem. I.e., the weight vector
#' \eqn{(w_1, w_2)} of an edge \eqn{(i, j)} is substituted with a weighted
#' sum \eqn{\lambda_i w_1 + (1 - \lambda_i) w_2} for different weights \eqn{\lambda_i \in [0, 1]}.
#'
#' @note Note that multi-criteria Prim may be trapped by intractability of the
#' multi-crtiteria spanning tree problem.
#'
#' @references J. D. Knowles and D. W. Corne, "A comparison of encodings and
#' algorithms for multiobjective minimum spanning tree problems," in Proceedings
#' of the 2001 Congress on Evolutionary Computation (IEEE Cat. No.01TH8546),
#' vol. 1, 2001, pp. 544–551 vol. 1.
#'
#' @param instance [\code{\link[grapherator]{grapherator}}]\cr
#'   Graph.
#' @return [\code{list}] List with component \code{pareto.front}.
#' @examples
#' g = genRandomMCGP(10)
#' res = mcMSTPrim(g)
#' print(res$pareto.front)
#' @family mcMST algorithms
#' @export
mcMSTPrim = function(instance) {
  if (inherits(instance, "grapherator"))
    instance = grapheratorToGraph(instance)
  n.weights = instance$getW()
  if (n.weights != 2L)
    stopf("mcMSTPrim: At the moment only bi-objective problems supported.")

  result = instance$doMCPrim()

  pareto.front = ecr::toParetoDf(result$weights)
  n.trees = nrow(pareto.front)

  # the trees are stored in a large matrix. Each two consecitive rows describe
  # one efficient tree
  pareto.set = lapply(seq(1, n.trees, by = 2L), function(i) {
    result$edges[i:(i + 1), , drop = FALSE]
  })

  res = list(
    pareto.set = pareto.set,
    pareto.front = pareto.front
  )
  return(ecr::filterDuplicated(res))
}
