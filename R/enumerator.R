#' @title Enumerate all solution candidates
#'
#' @description These functions enumerate all candidate solutions for
#' a certain combinatorial optimization problem, e.g., all permutations
#' for a TSP or all Pruefer-codes for a MST problem.
#'
#' @param n [\code{integer(1)}]\cr
#'   Instance size.
#' @return [\code{matrix}] Each row contains a candidate solution.
#' @example
#' sols = enumerateTSP(4L)
#' @rdname enumerators
#' @export
enumerateTSP = function(n) {
  n = asInt(n, lower = 3L)
  requirePackages("gtools", why = "rmoco::enumerateTSP")
  perms = gtools::permutations(n - 1L, n - 1L, v = 2:n)
  # drop duplicated permutations (symmetric TSP)
  perms = perms[1:(nrow(perms) / 2), , drop = FALSE]
  perms = cbind(matrix(1L, nrow = nrow(perms)), perms)
  return(perms)
}

#' @rdname enumerators
#' @export
enumerateMCMST = function(n) {
  n = asInt(n, lower = 3L)
  requirePackages("gtools", why = "rmoco::enumerateMCMST")
  gtools::permutations(n, n - 2L, seq_len(n), repeats.allowed = TRUE)
}
