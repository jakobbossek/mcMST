#' @title Subgraph EMOA for the multi-criteria MST problem.
#'
#' @description Evolutionary multi-objective algorithm to solve the
#' multi-objective minimum spanning tree problem. The algorithm relies
#' to mutation only to generate offspring. The package contains the subgraph mutator
#' (see \code{\link{mutSubgraphMST}}) or a simple one-edge exchange mutator
#' (see \code{\link{mutEdgeExchange}}). Of course, the user may use any
#' custom mutator which operators on edge lists as well
#' (see \code{\link[ecr]{makeMutator}}).
#'
#' @references Bossek, J., and Grimme, C. A Pareto-Beneficial Sub-Tree Mutation
#' for the Multi-Criteria Minimum Spanning Tree Problem. In Proceedings of the
#' 2017 IEEE Symposium Series on Computational Intelligence (2017). (accepted)
#'
#' @template arg_instance
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
#' @param ref.point [\code{numeric(n.objectives) | NULL}]\cr
#'   Reference point for hypervolume computation used for logging.
#'   If \code{NULL} the sum of the \eqn{n} largest edges in each objective
#'   is used where \eqn{n} is the number of nodes of \code{instance}.
#'   This is an upper bound for the size of each spanning tree
#'   with \eqn{(n-1)} edges.
#' @param max.iter [\code{integer(1)}]\cr
#'   Maximal number of iterations.
#'   Default is \code{100}.
#' @param ... [\code{any}]\cr
#'   Further parameters passed to mutator.
#' @template ret_ecrresult
#' @examples
#' inst = genRandomMCGP(10)
#' res = mcMSTEmoaBG(inst, mu = 20L, max.iter = 100L)
#' print(res$pareto.front)
#' print(tail(getStatistics(res$log)))
#' @family mcMST EMOAs
#' @family mcMST algorithms
#' @seealso Mutators \code{\link{mutSubgraphMST}} and \code{\link{mutEdgeExchange}}
#' @export
mcMSTEmoaBG = function(instance,
  mu, lambda = mu,
  mut = NULL,
  selMating = NULL, selSurvival = ecr::selNondom,
  ref.point = NULL,
  max.iter = 100L,
  ...) {

  # get number of nodes
  n = grapherator::getNumberOfNodes(instance)
  n.objectives = grapherator::getNumberOfWeights(instance)

  force(instance)

  # default is our subgraph mutator
  if (is.null(mut))
    mut = setup(mutSubgraphMST, instance = instance, ...)

  if (is.null(ref.point))
    ref.point = getReferencePoint(instance)

  fitness.fun = function(edgelist, instance) {
    getWeight(instance, edgelist)
  }

  # now generate an initial population, i.e.,
  # a list of edge lists
  population = lapply(1:mu, function(i) {
    #pcode = sample(1:n, n - 2L, replace = TRUE)
    #prueferToEdgeList(pcode)
    getRandomSpanningTree(instance)
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

getRandomSpanningTree = function(g) {
  adj.mat = grapherator::getAdjacencyMatrix(g)
  n = grapherator::getNumberOfNodes(g)

  # construct random distance matrix
  dmat = matrix(100 * runif(n * n), ncol = n, nrow = n)

  # stick to adjacency structure
  if (!is.null(adj.mat))
    dmat[!adj.mat] = 1e7
  else
    diag(dmat) = 1e7

  nodes = 1:n

  mstres = vegan::spantree(d = dmat)

  edge.list = matrix(
    c(nodes[2:n], nodes[mstres$kid]),
    byrow = TRUE, nrow = 2L)
  return(edge.list)
}
