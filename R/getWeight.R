#' Get the overall costs/weight of a suggraph given its edgelist.
#'
#' @param instance [\code{any}]\cr
#'   Multi-obective graph problem.
#' @template arg_edgelist
#' @return [numeric(2)] Weight vector.
#' @export
getWeight = function(instance, edgelist) {
  assertClass(instance, "mcGP")
  assertMatrix(edgelist)
  m = ncol(edgelist)

  # finally compute weights
  #FIXME: generalize to o > 2. Requires internal changes
  # in reprensentation of instances, i.e., list of matrizes
  w1 = w2 = 0
  for (i in seq_len(m)) {
    w1 = w1 + instance$w1[edgelist[1L, i], edgelist[2L, i]]
    w2 = w2 + instance$w2[edgelist[1L, i], edgelist[2L, i]]
  }

  return(c(w1, w2))
}
