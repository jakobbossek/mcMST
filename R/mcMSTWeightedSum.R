#' @title Multi-Objective Prim algorithm.
#'
#' @description Approximates the Pareto-optimal mcMST front of a multi-objective
#' graph problem by iteratively applying Prim's algorithm for the single-objective
#' MST problem to a scalarized version of the problem. I.e., the weight vector
#' \eqn{(w_1, w_2)} of an edge \eqn{(i, j)} is substituted with a weighted
#' sum \eqn{\lambda_i w_1 + (1 - \lambda_i) w_2} for different weights \eqn{\lambda_i \in [0, 1]}.
#'
#' @note Note that this procedure can only find socalled supported efficient
#' solutions, i.e., solutions on the convex hull of the Pareto-optimal front.
#'
#' @references J. D. Knowles and D. W. Corne, "A comparison of encodings and
#' algorithms for multiobjective minimum spanning tree problems," in Proceedings
#' of the 2001 Congress on Evolutionary Computation (IEEE Cat. No.01TH8546),
#' vol. 1, 2001, pp. 544–551 vol. 1.
#'
#' @param instance [\code{\link[grapherator]{grapherator}}]\cr
#'   Graph.
#' @param n.lambdas [\code{integer(1) | NULL}]\cr
#'   Number of weights to generate. The weights are generated equdistantly
#'   in the interval \eqn{[0, 1]}.
#' @param lambdas [\code{numerci}]\cr
#'   Vector of weights. This is an alternative to \code{n.lambdas}.
#' @return [\code{list}] List with component \code{pareto.front}.
#' @examples
#' g = genRandomMCGP(30)
#' res = mcMSTWeightedSum(g, n.lambdas = 50)
#' print(res$pareto.front)
#' @family mcMST algorithms
#' @export
mcMSTWeightedSum = function(instance, n.lambdas = NULL, lambdas = NULL, weightSamplingFun = defaultWeightSampling) {
  if (inherits(instance, "grapherator"))
    instance = grapheratorToGraph(instance)
  n.weights = instance$getW()
  if (n.weights < 2L)
    BBmisc::stopf("[mcMST::mcMSTWeightedSum] Requires at least two weights per edge.")
  if (is.null(n.lambdas) & is.null(lambdas))
    BBmisc::stopf("[mcMST::mcMSTWeightedSum] At least n.lambdas or lambdas must be set.")
  if (is.null(lambdas) & !is.null(n.lambdas)) {
    n.lambdas = asInt(n.lambdas)
    lambdas = seq(0, 1, length.out = n.lambdas)
  }
  assertNumeric(lambdas, any.missing = FALSE, all.missing = FALSE)

  pareto.set = vector(mode = "list", length = length(lambdas))
  pareto.front = matrix(0, ncol = length(lambdas), nrow = n.weights)

  n.lambdas = length(lambdas)
  # now apply scalarization and compute supported efficient solution(s)
  for (k in seq_len(n.lambdas)) {
    tree = instance$getMSTByWeightedSumScalarization(lambdas[k])
    pareto.front[, k] = tree$getSumOfEdgeWeights()
    pareto.set[[k]] = tree$toEdgeList()
  }

  # filter dominated solutions
  idx.nondom = ecr::which.nondominated(pareto.front)
  pareto.front = pareto.front[, idx.nondom, drop = FALSE]
  pareto.set = pareto.set[idx.nondom]

  res = list(
    pareto.set = pareto.set,
    pareto.front = ecr::toParetoDf(pareto.front)
  )
  return(ecr::filterDuplicated(res))
}

