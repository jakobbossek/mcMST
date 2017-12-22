#' @title Enumerate all Pareto-optimal solutions.
#'
#' @description Function which expects a problem instance of a combinatorial optimization
#' problem (e.g., MST), a multi-objective function and a solution enumerator, i.e., a function
#' which enumerates all possible solutions (e.g., all Pruefer codes in case of a
#' MST problem) and  determines both the Pareto front and Pareto set by
#' exhaustive enumeration.
#'
#' @note This method exhaustively enumerates all possible solutions
#' of a given multi-objective combinatorial optimization problem. Thus,
#' it is limited to small input size due to combinatorial explosion.
#'
#' @param instance [any]\cr
#'   Problem instance.
#' @param obj.fun [\code{function(solution, instance)}]\cr
#'   Objective function which expects a numeric vector \code{solution} encoding a
#'   solution candidate and a problem instance \code{instance}. The function should
#'   return a numeric vector of length \code{n.objectives}.
#' @param enumerator.fun [\code{function(n)}]\cr
#'   Function to exhaustively generate all possible candidate solutions.
#'   Expects a single integer value n, i.e., the instance size, e.g., the
#'   number of nodes for a graph problem.
#' @param n.objectives [\code{integer(1)}]\cr
#'   Number of objectives of problem.
#' @param simplify [\code{logical(1)}]\cr
#'   Should pareto set be simplified to matrix?
#'   This will only be done if all elements are of the same length. Otherwise
#'   the parameter will be ignored.
#'   Default is \code{TRUE}.
#' @return [\code{list}] List with elements \code{pareto.set} (matrix of Pareto-optimal solutions)
#' and \code{pareto.front} (matrix of corresponding weight vectors).
#' @examples
#' # here we enumerate all Pareto-optimal solutions of a bi-objective mcMST problem
#' # we use the Pruefer-code enumerator. Thus, we need to define an objective
#' # function, which is able to handle this type of endcoding
#' objfunMCMST = function(pcode, instance) {
#'   getWeight(instance, prueferToEdgeList(pcode))
#' }
#'
#' # next we generate a random bi-objective graph
#' g = genRandomMCGP(5L)
#'
#' # ... and finally compute the exact front of g
#' res = getExactFront(g, obj.fun = objfunMCMST, enumerator.fun = enumerateMST, n.objectives = 2L)
#' \dontrun{
#' plot(res$pareto.front)
#' }
#' @export
getExactFront = function(instance, obj.fun, enumerator.fun, n.objectives, simplify = TRUE) {
  assertClass(instance, "grapherator")
  assertFunction(obj.fun)
  assertFunction(enumerator.fun, args = "n")
  n.objectives = asInt(n.objectives, lower = 2L)
  assertFlag(simplify)

  n.nodes = grapherator::getNumberOfNodes(instance)
  if (n.nodes > 10L)
    warningf("Doh! This may take some time.")

  # allocate really large vector of permutations
  pp = enumerator.fun(n.nodes)
  n.sols = if (is.matrix(pp)) nrow(pp) else length(pp)

  # convert to list in any case
  matrix.pp = is.matrix(pp)
  if (matrix.pp)
    pp = lapply(1:nrow(pp), function(i) pp[i, ])

  pareto.set = list()
  pareto.front = matrix(nrow = n.objectives, ncol = 0L)

  n.step = 1000L

  # now perform kind of a sliding window approach. This is done to avoid doing
  # non-domination check for n^(n-2) elements.
  i = 1L
  while (i < n.sols) {
    # which elements to check next
    j = min(i + n.step, n.sols)

    # compute objectives
    weights = lapply(pp[i:j], obj.fun, instance)
    weights = do.call(cbind, weights)

    # update sets
    pareto.set = c(pareto.set, pp[i:j])
    pareto.front = cbind(pareto.front, weights)

    # filter dominated
    idx.nondom = ecr::which.nondominated(pareto.front)
    pareto.set = pareto.set[idx.nondom]
    pareto.front = pareto.front[, idx.nondom, drop = FALSE]
    i = j + 1L
  }

  if (simplify & matrix.pp)
    pareto.set = do.call(rbind, pareto.set)

  return(list(
    pareto.set = pareto.set,
    pareto.front = pareto.front)
  )
}

getExactFrontMCMST = function(instance, ...) {
  objfunMCMST = function(pcode, instance) {
    #print(pcode)
    getWeight(instance, prueferToEdgeList(pcode))
  }
  getExactFront(instance, obj.fun = objfunMCMST,
    enumerator.fun = enumerateMST, n.objectives = 2L, ...)
}
