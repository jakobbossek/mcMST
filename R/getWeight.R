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

  # finally compute weights
  weights1 = graph$weights[[1L]]
  weights2 = graph$weights[[2L]]
  w1 = w2 = 0
  for (i in seq_len(m)) {
    w1 = w1 + weights1[edgelist[1L, i], edgelist[2L, i]]
    w2 = w2 + weights2[edgelist[1L, i], edgelist[2L, i]]
  }

  return(c(w1, w2))
}
