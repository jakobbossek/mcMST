# Determine exact Pareto-front.
#
# Note: the instance size needs to be really low,
# since we exhaustively enumerate all solutions.
#
# @param instance [any]
#   Problem instance.
# @param obj.fun [function]
#   Objective function which expects a numeric vector and an instance.
# @param enumerator.fun [function(n)]
#   Function to exhaustively generate all possible solutions.
#   Expects a single integer value n, i.e., the instance size, e.g., the
#   number of nodes for a graph problem.
# @return [list] List with elements pareto.set (matrix of Pruefer codes) and
# pareto.front (matrix of weight vectors).
getExactFront = function(instance, obj.fun, enumerator.fun, n.objectives) {
  assertFunction(obj.fun)
  assertFunction(enumerator.fun, args = "n")
  n.objectives = asInt(n.objectives, lower = 2L)


  n = instance$n
  if (n > 10L)
    warningf("Doh! This may take some time.")

  # allocate really large vector of permutations
  pp = enumerator.fun(n)
  n.sols = nrow(pp)
  len.sol = ncol(pp)

  pareto.set = matrix(ncol = len.sol, nrow = 0L)
  pareto.front = matrix(nrow = n.objectives, ncol = 0L)

  n.step = 1000L

  # now perform kind of a sliding window approach. This is done to avoid doing
  # non-domination check for n^(n-2) elements.
  i = 1L
  while (i < n.sols) {
    j = min(i + n.step, n.sols)
    weights = apply(pp[i:j, , drop = FALSE], 1L, obj.fun, instance)

    pareto.set = rbind(pareto.set, pp[i:j, , drop = FALSE])
    pareto.front = cbind(pareto.front, weights)

    idx.nondom = ecr::which.nondominated(pareto.front)
    pareto.set = pareto.set[idx.nondom, , drop = FALSE]
    pareto.front = pareto.front[, idx.nondom, drop = FALSE]
    i = j + 1L
  }

  return(list(pareto.set = pareto.set, pareto.front = pareto.front))
}
