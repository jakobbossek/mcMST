#' @title Sample weights
#'
#' @description Sample random weights \eqn{\lambda_1, ... \lambda_n, \sum_{i=1}^{n} \lambda_i = 1}
#' for weighted-sum scalarization.
#'
#' @param n [\code{integer(1)}]\cr
#'   Number of weights to sample.
#' @return [\code{numeric}] Weight vector.
#' @examples
#' sampleWeights(2)
#'
#' weights = replicate(10, sampleWeights(3L))
#' colSums(weights)
#' @export
sampleWeights = function(n) {
  n = asInt(n, lower =  2L)
  lambdas = numeric(n)
  max.val = 1
  for (i in 1:(n - 1L)) {
    lambdas[i] = runif(1L, max = max.val)
    max.val = max.val - lambdas[i]
  }
  lambdas[n] = max.val
  return(lambdas)
}

#' @title Scalarize weight matrizes.
#'
#' @description Given a list of weight matrizes \code{weight.mats} and a vector
#' of numeric weights, the function returns a single weight matrix. Each component
#' of the resulting matrix is the weighted sum of the corresponding components of
#' the weight matrizes passed.
#'
#' @param weight.mats [\code{list}]\cr
#'   List of weight matrizes.
#' @param lambdas [\code{numeric}]\cr
#'   Vector of weights.
#' @return [\code{matrix}]
#' @export
scalarizeWeights = function(weight.mats, lambdas) {
  assertList(weight.mats, types = "matrix")
  assertNumeric(lambdas, lower = 0, upper = 1, any.missing = FALSE, all.missing = FALSE)
  n = length(weight.mats)
  if (n != length(lambdas))
    stopf("scalarize: number of objectives not equal to number of lambda values!")
  # accumulated values
  acc = matrix(0, nrow = nrow(weight.mats[[1L]]), ncol = ncol(weight.mats[[1L]]))
  for (i in 1:n) {
    acc = acc + lambdas[i] * weight.mats[[i]]
  }
  return(acc)
}
