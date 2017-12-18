#' @title Compute similarity matrix.
#'
#' @description Given a list of objects and a function which computes a similarity
#' measure between two objects of the list, \code{computeSimilarity} returns a
#' similarity matrix.
#'
#' @param set [\code{list}]\cr
#'   List of objects.
#' @param sim.fun [\code{function(x, y, ...)}]\cr
#'   Function which expects two objects \code{x} and \code{y} as first and second
#'   arguments and returns a scalar value.
#' @param ... [any]\cr
#'   Passed down to \code{sim.fun}.
#' @return [\code{matrix(n, n)}] \eqn{(n,n)} matrix with \eqn{n} being the length of \code{set}.
#' @export
computeSimilarityMatrix = function(set, sim.fun, ...) {
  assertList(set)
  assertFunction(sim.fun)

  sim.fun = Vectorize(sim.fun, vectorize.args = c("x", "y"), SIMPLIFY = TRUE)

  sim.mat = outer(set, set, sim.fun, ...)
  rownames(sim.mat) = colnames(sim.mat) = paste0("S", 1:nrow(sim.mat))
  return(sim.mat)
}
