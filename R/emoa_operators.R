#' @title Uniform mutation for Pruefer code representation.
#'
#' @description \code{mutUniformPruefer} replaces each component of a Pruefer code of length n - 2
#' with probability \code{p} with a random node number between 1 and n.
#'
#' @param ind [\code{integer}]\cr
#'   Pruefer code.
#' @param p [\code{numeric(1)}]\cr
#'   Probability of mutation of each component of \code{ind}.
#'   Default is \code{1 / length(ind)}.
#' @return [\code{integer}] Mutated Pruefer code.
#' @family mcMST EMOA mutators
#' @seealso Evolutionary multi-objective algorithm \code{\link{mcMSTEmoaZhou}}
#' @export
mutUniformPruefer = makeMutator(
  mutator = function(ind, p = 1 / length(ind)) {
    n = length(ind)
    nrs = 1:(n + 2L)
    r = runif(n) < p
    if (sum(r) > 0)
      ind[r] = sample(nrs, sum(r), replace = TRUE)
    return(ind)
  },
  supported = "custom")

oneEdgeExchange = function(edgelist, edge.id, instance) {
  n = grapherator::getNumberOfNodes(instance)

  the.edge.idx = edge.id
  # get incident nodes
  node1 = edgelist[1L, the.edge.idx]
  node2 = edgelist[2L, the.edge.idx]

  getRandom = function(set) {
    if (length(set) == 1L)
      return(set)
    return(sample(set, 1L))
  }

  # build igraph from edgelist
  tmp = igraph::graph_from_edgelist(t(edgelist), directed = FALSE)

  # now delete selected edges
  tmp = igraph::delete_edges(tmp, collapse(edgelist[, the.edge.idx], sep = "|"))

  # ... and search for components
  components = igraph::components(tmp, mode = "weak")
  comp.node1 = components$membership[node1]
  comp.node2 = components$membership[node2]

  nodes1 = which(components$membership == comp.node1)
  nodes2 = which(components$membership == comp.node2)

  getRandomConnectingEdge = function(nodes1, nodes2) {
    adj.mat = grapherator::getAdjacencyMatrix(instance)

    # reduce to part of adjacency matrix
    adj.mat2 = adj.mat[nodes1, nodes2, drop = FALSE]

    # get matrix rows and cols
    pos.edges = which(adj.mat2, arr.ind = TRUE)
    if (!is.matrix(pos.edges))
      pos.edges = matrix(pos.edges, nrow = 1L)

    # sample an edge at random
    idx.newnodes = pos.edges[getRandom(1:nrow(pos.edges)), ]

    return(c(nodes1[idx.newnodes[1L]], nodes2[idx.newnodes[2L]]))
  }

  # replace edge with new edge
  edgelist[, the.edge.idx] = getRandomConnectingEdge(nodes1, nodes2)
  #checkValidity(edgelist, instance)

  return(edgelist)
}

# Replace each edge with a given probability
edgeExchange = function(edgelist, p = 1 / ncol(edgelist), instance = NULL) {
  # get number of tree edges
  m = ncol(edgelist)

  # replace each edge with prob p
  for (edge.id in 1:m) {
    if (runif(1L) < p) {
      # sample random edge in tree
      edgelist = oneEdgeExchange(edgelist, edge.id, instance = instance)
    }
  }
  return(edgelist)
}

#' @title One-edge-exchange mutator for edge list representation of spanning trees.
#'
#' @description Each edge is replaced with another feasible edge with probability p.
#' By default p = 1/m where m is the number of edges, i.e., in expectation one edge
#' is replaced. The operators maintains the spanning tree property, i.e., the resulting
#' edge list is indeed the edge list of a spanning tree.
#'
#' @param ind [\code{matrix(2, m)}]\cr
#'   Matrix of edges (each column is one edge).
#' @param p [\code{numeric(1)}]\cr
#'   Probability of edge exchange.
#'   Default is \code{1 / ncol(ind)}.
#' @template arg_instance
#' @return [\code{matrix(2, m)}] Mutated edge list.
#' @family mcMST EMOA mutators
#' @seealso Evolutionary multi-objective algorithm \code{\link{mcMSTEmoaBG}}
#' @export
mutEdgeExchange = makeMutator(
  mutator = function(ind, p = 1 / ncol(ind), instance = NULL) {
    edgeExchange(ind, p, instance)
  },
  supported = "custom")


#' @title k-edge-exchange mutator for edge list representation of spanning trees.
#'
#' @description Let \eqn{m} be the number of spanning tree edges. Then, the operator
#' selects \eqn{1 \leq k \leq m} edges randomly and replaces each of the \eqn{k}
#' edges with another feasible edge.
#'
#' @inheritParams mutEdgeExchange
#' @param k [\code{integer(1)}]\cr
#'   Number of edges to swap.
#' @family mcMST EMOA mutators
#' @seealso Evolutionary multi-objective algorithm \code{\link{mcMSTEmoaBG}}
#' @export
mutKEdgeExchange = makeMutator(
  mutator = function(ind, k = 1L, drop.largest = FALSE, instance = NULL) {
    instance$getMSTByEdgeExchange(ind, k, drop.largest)
  },
  supported = "custom")

#' @title Subgraph-mutator for edge list representation.
#'
#' @description \code{mutSubgraphMST} selects a random edge e = (u, v) and traverses
#' the tree starting form u and v respectively until a connected subtree of at most
#' \code{sigma} edges is selected. Then the subtree is replaced with the optimal spanning subtree
#' regarding one of the objectives with equal probability.
#'
#' @param ind [\code{matrix(2, m)}]\cr
#'   Matrix of edges (each column is one edge).
#' @param sigma [\code{integer()}]\cr
#'   Upper bound for the size of the selected subtree.
#' @param scalarize [\code{logical(1)}]\cr
#'   Should a scalarized version of the the subproblem be solved?
#'   If \code{TRUE}, a random weight \eqn{\lambda \in [0,1]} is sampled
#'   from a \code{U[0, 1]}-distribution. Next, a weighted sum
#'   scalarization \eqn{\lambda \cdot c_1 + (1 - \lambda) \cdot c_2}
#'   of the subproblem is solved.
#'   Default is \code{FALSE}, i.e., the single-objective subproblem is
#'   solved. The objective to optimize for is sampled with equal probability.
#' @template arg_instance
#' @return [\code{matrix(2, m)}] Mutated edge list.
#' @family mcMST EMOA mutators
#' @seealso Evolutionary multi-objective algorithm \code{\link{mcMSTEmoaBG}}
#' @export
mutSubgraphMST = makeMutator(
  mutator = function(ind, sigma = floor(ind$getV() / 2L), scalarize = FALSE, altWeightGen = FALSE, instance = NULL) {
    n.select = sample(3:sigma, 1L)
    instance$getMSTBySubgraphMutation(ind, n.select, scalarize = scalarize, altWeightGen = altWeightGen)
  },
  supported = "custom"
)

#' @title Forest-mutator for edge list representation.
#'
#' @description \code{mutForestMST} drops k edges randomly. In consequence the
#' tree is decomposed into k+1 subtrees (forest). Now the operator reconnects the
#' subtrees by constructing a minimum spanning tree between the components.
#'
#' @param ind [\code{matrix(2, m)}]\cr
#'   Matrix of edges (each column is one edge).
#' @param sigma [\code{integer()}]\cr
#'   Upper bound for number of dropped edges.
#' @param scalarize [\code{logical(1)}]\cr
#'   Should a scalarized version of the the subproblem be solved?
#'   If \code{TRUE}, a random weight \eqn{\lambda \in [0,1]} is sampled
#'   from a \code{U[0, 1]}-distribution. Next, a weighted sum
#'   scalarization \eqn{\lambda \cdot c_1 + (1 - \lambda) \cdot c_2}
#'   of the subproblem is solved.
#'   Default is \code{FALSE}, i.e., the single-objective subproblem is
#'   solved. The objective to optimize for is sampled with equal probability.
#' @template arg_instance
#' @return [\code{matrix(2, m)}] Mutated edge list.
#' @family mcMST EMOA mutators
#' @seealso Evolutionary multi-objective algorithm \code{\link{mcMSTEmoaBG}}
#' @export
mutSubforestMST = makeMutator(
  mutator = function(ind, sigma = floor(ind$getV() / 2), scalarize = FALSE, instance = NULL) {
    n.select = sample(1:sigma, 1L)
    instance$getMSTBySubforestMutation(ind, n.select)
  },
  supported = "custom"
)
