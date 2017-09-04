#' @return [\code{\link[ecr]{ecr_result}}] List of type \code{\link[ecr]{ecr_result}}
#'  with the following components:
#'  \describe{
#'   \item{task}{The \code{ecr_optimization_task}.}
#'   \item{log}{Logger object.}
#'   \item{pareto.idx}{Indizes of the non-dominated solutions in the last population.}
#'   \item{pareto.front}{(n x d) matrix of the approximated non-dominated front where n
#'   is the number of non-dominated points and d is the number of objectives.}
#'   \item{pareto.set}{Matrix of decision space values resulting with objective values
#'   given in pareto.front.}
#'   \item{last.population}{Last population.}
#'   \item{message}{Character string describing the reason of termination.}
#'  }
