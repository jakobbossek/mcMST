#' @title Subgraph EMOA for the multi-criteria MST problem.
#'
#' @description Evolutionary multi-objective algorithm to solve the
#' multi-objective minimum spanning tree problem. The algorithm relies
#' to mutation only to generate offspring employing the subgraph mutator
#' (see \code{\link{mutSubgraphMST}}).
#'
#' @note This algorithm performs minimization in all objectives.
#'
#' @param instance [\code{any}]\cr
#'   Multi-objective graph problem.
#' @param n.objectives [\code{integer(1)}]\cr
#'   Number of objectives.
#' @param mu [\code{integer(1)}]\cr
#'   Population size.
#' @param lambda [\code{integer(1)}]\cr
#'   Number of offspring generated in each generation.
#'   Default is \code{mu}.
#' @param mut [\code{ecr_mutator}]\cr
#'   Mutation operator.
#'   Default is \code{\link{mutSubgraphMST}}.
#' @param selMating [\code{ecr_selector}]\cr
#'   Mating selector.
#'   Default is \code{\link[ecr]{selSimple}}.
#' @param selSurvival [\code{ecr_selector}]\cr
#'   Survival selector.
#'   Default is \code{link[ecr]{selNondom}}.
#' @param ref.point [\code{numeric(n.objectives)}]\cr
#'   Reference point for Hypervolume computation.
#' @param max.iter [\code{integer(1)}]\cr
#'   Maximal number of iterations.
#'   Default is \code{100}.
#' @return [\code{\link[ecr]{ecr_result}}]
#' @examples
#' inst = genRandomMCGP(10)
#' res = mcMSTEmoaBG(inst, mu = 20L, ref.point = c(1e5, 1e5), max.iter = 100L)
#' print(res$pareto.front)
#' print(tail(getStatistics(res$log)))
#' @family mcMST EMOAs
#' @export
mcMSTEmoaBG = function(instance, n.objectives = 2L,
  mu, lambda = mu,
  mut = NULL,
  selMating = NULL, selSurvival = ecr::selNondom,
  ref.point,
  max.iter = 100L) {

  # get number of nodes
  n = instance$n.nodes

  force(instance)

  # default is our subgraph mutator
  if (is.null(mut))
    mut = setup(mutSubgraphMST, instance = instance)

  fitness.fun = function(edgelist, instance) {
    getWeight(instance, edgelist)
  }

  # now generate an initial population, i.e.,
  # a list of edge lists
  population = lapply(1:mu, function(i) {
    pcode = sample(1:n, n - 2L, replace = TRUE)
    prueferToEdgelist(pcode)
  })

  res = ecr(fitness.fun = fitness.fun, n.objectives = n.objectives,
    mu = mu, lambda = lambda, survival.strategy = "plus", representation = "custom",
    initial.solutions = population,
    survival.selector = selSurvival, parent.selector = selMating,
    mutator = mut, p.mut = 1,
    log.stats = list(fitness = list("HV" = list(
      fun = computeHV,
      pars = list(ref.point = ref.point)))),
    terminators = list(stopOnIters(max.iter)),
    instance = instance)

  return(res)
}
