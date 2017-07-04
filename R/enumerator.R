enumerateTSP = function(n) {
  n = asInt(n, lower = 3L)
  requirePackages("gtools", why = "enumerateTSP")
  perms = gtools::permutations(n - 1L, n - 1L, v = 2:n)
  # drop duplicated permutations (symmetric TSP)
  perms = perms[1:(nrow(perms) / 2), , drop = FALSE]
  perms = cbind(matrix(1L, nrow = nrow(perms)), perms)
  return(perms)
}

enumerateMCMST = function(n) {
  n = asInt(n, lower = 3L)
  requirePackages("gtools", why = "enumerateMCMST")
  gtools::permutations(n, n - 2L, seq_len(n), repeats.allowed = TRUE)
}
