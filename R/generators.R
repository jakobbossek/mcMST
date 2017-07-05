#' @title Generate a random bi-criteria graph instance.
#'
#' @description No topology is defined. The instance is composed of two
#' symmetric weight matrices.
#'
#' @param n [\code{integer(1)}]\cr
#'   Instance size, i.e., number of nodes.
#' @return [\code{mcGP}]
#' @examples
#' inst = genRandomMCGP(6L)
#' @export
genRandomMCGP = function(n) {
  n = asInt(n, lower = 2L)
  w1 = matrix(runif(n^2, min = 10, max = 100), nrow = n)
  w1[lower.tri(w1)] = t(w1)[lower.tri(t(w1))]
  w2 = matrix(runif(n^2, min = 10, max = 50), nrow = n)
  w2[lower.tri(w2)] = t(w2)[lower.tri(t(w2))]
  diag(w1) = diag(w2) = 0
  BBmisc::makeS3Obj(classes = "mcGP", n = n, w1 = w1, w2 = w2)
}
