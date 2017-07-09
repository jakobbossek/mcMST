#' Get the overall costs/weight of a subgraph given its edgelist.
#'
#' @template arg_mcGP
#' @template arg_edgelist
#' @return [\code{numeric(2)}] Weight vector.
#' @export
getWeight = function(graph, edgelist) {
  assertClass(graph, "mcGP")
  assertMatrix(edgelist)
  m = ncol(edgelist)

  n.weights = graph$n.weights

  # finally compute weights
  ws = numeric(n.weights)
  #FIXME: inefficient
  for (i in seq_len(m)) {
    for (j in seq_len(n.weights)) {
      ws[j] = ws[j] + graph$weights[[j]][edgelist[1L, i], edgelist[2L, i]]
    }
  }

  return(ws)
}
