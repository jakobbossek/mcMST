#' Get the overall costs/weight of a subgraph given its edgelist.
#'
#' @template arg_mcGP
#' @template arg_edgelist
#' @param obj.type [\code{character(1)}]\cr
#'   How to aggregate edge weights?
#'   Possible values are \dQuote{sum} for sum objectives and \dQuote{bottleneck}
#'   for bottleneck/min-max objectives.
#'   Default is \dQuote{sum}.
#' @return [\code{numeric(2)}] Weight vector.
#' @examples
#' # generate a random bi-objective graph
#' g = genRandomMCGP(5)
#'
#' # generate a random Pruefer code, i.e., a random spanning tree of g
#' pcode = sample(1:5, 3, replace = TRUE)
#'
#' getWeight(g, prueferToEdgeList(pcode))
#' getWeight(g, prueferToEdgeList(pcode), obj.type = "bottleneck")
#' @export
getWeight = function(graph, edgelist, obj.type = "sum") {
  assertClass(graph, "mcGP")
  assertMatrix(edgelist)
  assertChoice(obj.type, choices = c("sum", "bottleneck"))
  m = ncol(edgelist)

  n.weights = graph$n.weights

  # get edge weights (one column for each weight)
  edge.weights = getWeights(graph, edgelist)
  aggr.fun = ifelse(obj.type == "sum", sum, max)
  obj.vec = apply(edge.weights, 1L, aggr.fun)
  return(obj.vec)
}


getWeights = function(graph, edgelist) {
  assertClass(graph, "mcGP")
  assertMatrix(edgelist)
  m = ncol(edgelist)

  n.weights = graph$n.weights

  # finally compute weights
  ws = matrix(NA, ncol = m, nrow = n.weights)

  #FIXME: inefficient
  for (i in seq_len(m)) {
    for (j in seq_len(n.weights)) {
      ws[j, i] = graph$weights[[j]][edgelist[1L, i], edgelist[2L, i]]
    }
  }
  return(ws)
}
