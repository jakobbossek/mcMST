#' @title Generate a bi-criteria graph with uniformly randomly distribted edge weights.
#'
#' @description No topology is defined. The instance is composed of two
#' symmetric weight matrices. The first weight is drawn independently at
#' random from a \eqn{\mathcal{R}[10, 100]} distribution, the second one
#' from a \eqn{\mathcal{R}[10, 50]} distribution (see references).
#'
#' @note This is a simple wrapper around the much more flexible graph generation
#' system (see, e.g., \code{\link{mcGP}}).
#'
#' @references
#' Zhou, G. and Gen, M. Genetic Algorithm Approach on Multi-Criteria
#' Minimum Spanning Tree Problem. In: European Journal of Operational Research (1999).
#'
#' Knowles, JD & Corne, DW 2001, A comparison of encodings and algorithms for multiobjective
#' minimum spanning tree problems. in Proceedings of the IEEE Conference on Evolutionary
#' Computation, ICEC|Proc IEEE Conf Evol Comput Proc ICEC. vol. 1, Institute of Electrical
#' and Electronics Engineers , pp. 544-551, Congress on Evolutionary Computation 2001,
#' Soul, 1 July.
#'
#' @param n [\code{integer(1)}]\cr
#'   Instance size, i.e., number of nodes.
#' @return [\code{mcGP}]
#' @examples
#' g = genRandomMCGP(10L)
#' \dontrun{
#' pl = plot(g)
#' }
#' @export
genRandomMCGP = function(n) {
  n = asInt(n, lower = 2L)
  g = mcGP(lower = 10, upper = 100)
  g = addCoordinates(g, n = n, generator = coordGrid)
  g = addWeights(g, method = "random", weight.fun = runif, symmetric = TRUE, min = 10, max = 100)
  g = addWeights(g, method = "random", weight.fun = runif, symmetric = TRUE, min = 10, max = 50)
  return(g)
}
