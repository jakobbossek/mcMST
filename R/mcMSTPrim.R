#' @title Multi-Objective Prim algorithm.
#'
#' @description Approximates the Pareto-optimal mcMST front of a multi-objective
#' graph problem by iteratively applying Prim's algorithm for the single-objective
#' MST problem to a scalarized version of the problem. I.e., the weight vector
#' \eqn{(w_1, w_2)} of an edge \eqn{(i, j)} is substituted with a weighted
#' sum \eqn{\lambda_i w_1 + (1 - \lambda_i) w_2} with weight \eqn{\lambda_i \in [0, 1]}
#' for different weights.
#'
#' @note Note that this procedure can only find socalled supported efficient
#' solutions, i.e., solutions on the convex hull of the Pareto-optimal front.
#'
#' @param instance [\code{mcGP}]\cr
#'   Multi-objective graph problem.
#' @param n.lambdas [\code{integer(1) | NULL}]\cr
#'   Number of weights to generate. The weights are generated equdistantly
#'   in the interval \eqn{[0, 1]}.
#' @param lambdas [\code{numerci}]\cr
#'   Vector of weights. This is an alternative to \code{n.lambdas}.
#' @return [\code{list}] List with component \code{pareto.front}.
#' g = genRandomMCGP(30)
#' res = mcMSTScalar(30, n.lambdas = 50)
#' print(res$pareto.front)
#' @export
#FIXME: generalize to > 2 objectives
mcMSTPrim = function(instance, n.lambdas = NULL, lambdas = NULL) {
  assertClass(instance, "mcGP")
  if (is.null(n.lambdas) & is.null(lambdas))
    stopf("mcMSTScalar: At least n.lambdas or lambdas must be set.")
  if (is.null(lambdas) & !is.null(n.lambdas)) {
    n.lambdas = asInt(n.lambdas)
    lambdas = seq(0, 1, length.out = n.lambdas)
  }
  assertNumeric(lambdas, any.missing = FALSE, all.missing = FALSE)

  #FIXME: also return pareto.set as Pruefer number (needs
  # transformation edgelistToPrueferCode)
  #pareto.set = matrix(0, nrow = )
  pareto.front = matrix(0, ncol = length(lambdas), nrow = instance$n.weights)

  # Helper function to build the weighted sum
  # of edge weights
  scalarize = function(instance, lambda) {
    weights = instance$weights
    assertList(weights, types = "matrix")
    assertNumber(lambda, lower = 0, upper = 1)

    lambda * weights[[1L]] + (1 - lambda) * weights[[2L]]
  }

  n = instance$n.nodes

  # now apply scalarization and compute supported efficient solution(s)
  for (k in 1:length(lambdas)) {
    weight.mat = scalarize(instance, lambda = lambdas[k])
    mst.res = vegan::spantree(d = weight.mat)
    nodes1 = 2:n
    nodes2 = mst.res$kid

    mst.costs = c(0, 0)
    for (i in 1:(n - 1)) {
      mst.costs = mst.costs + c(instance$weights[[1L]][nodes1[i], nodes2[i]],
        instance$weights[[2L]][nodes1[i], nodes2[i]])
    }
    pareto.front[, k] = mst.costs
  }
  return(list(pareto.front = pareto.front))
}
