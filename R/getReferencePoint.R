# @title Helper to approximate reference point for EMOA runs.
#
# @description Here we take the sum of the n largest edges (this is always an upper bound
# for each tree with (n-1) edges) in each objective.
# @param instance [\code{mcMST}]\cr
#   Problem instance
# @return [\code{numeric(2)}]
getReferencePoint = function(instance) {
  assertClass(instance, "mcGP")
  w1.sorted = sort(as.numeric(instance$weights[[1L]]), decreasing = TRUE)
  w2.sorted = sort(as.numeric(instance$weights[[2L]]), decreasing = TRUE)
  n = instance$n.nodes
  c(sum(w1.sorted[1:n]), sum(w2.sorted[1:n]))
}
