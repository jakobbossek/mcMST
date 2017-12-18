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

  #catf("End nodes: %i, %i", node1, node2)

  getRandom = function(set) {
    if (length(set) == 1L)
      return(set)
    return(sample(set, 1L))
  }

  #SLOW AS HELL! :)
  # getReachableNodes = function(edgelist, node) {
  #   rnodes = node
  #   nsel = 0
  #   nsel2 = 1
  #   while (nsel2 > nsel) {
  #     nsel = nsel2
  #     ridx = which(apply(edgelist, 2L, function(x) any(x %in% rnodes)))
  #     rnodes = c(unique(as.integer(edgelist[, ridx])), node)
  #     nsel2 = length(rnodes)
  #   }
  #   return(rnodes)
  # }

  #SLOW AS HELL! :)
  # getReachableNodes = function(edgelist, node) {
  #   nodes = node
  #   nodes2 = nodes
  #   while (length(nodes2) > 0L) {
  #     the.node = nodes2[1L]
  #     # check which edges are incident to current node
  #     is.adjacent = apply(edgelist, 2L, function(x) any(x == the.node))
  #     # get all unique nodes
  #     adjacent.nodes = unique(as.integer(edgelist[, is.adjacent]))
  #     nodes2 = setdiff(c(nodes2, adjacent.nodes), the.node)
  #     nodes = unique(c(nodes, adjacent.nodes))
  #     edgelist = edgelist[, !is.adjacent, drop = FALSE]
  #   }
  #   return(nodes)
  # }

  getReachableNodes = function(edgelist, node) {
    #print(instance)
    n = instance$n.nodes

    adj.list = replicate(n, c())
    for (i in 1:ncol(edgelist)) {
      n1 = edgelist[1L, i]
      n2 = edgelist[2L, i]
      #catf("n1: %i, n2: %i, n: %i", n1, n2, n)
      adj.list[[n1]] = c(adj.list[[n1]], n2)
      adj.list[[n2]] = c(adj.list[[n2]], n1)
    }
    adj.list.len = sapply(adj.list, length)

    #print(adj.list)

    queue = node
    qn = 1L
    visited = rep(FALSE, n)
    while (qn > 0L) {
      cur.node = queue[1L]
      #catf("Curnode: %i", cur.node)
      queue = queue[-1L]
      qn = qn - 1L
      if (visited[cur.node])
        next
      visited[cur.node] = TRUE
      queue = c(queue, adj.list[[cur.node]])
      qn = qn + adj.list.len[cur.node]
    }
    return(which(visited))
  }

  # now for each end node get all reachable nodes via DFS
  #print(edgelist)
  #print(edgelist)
  nodes1 = getReachableNodes(edgelist[, -the.edge.idx], node1)
  nodes2 = getReachableNodes(edgelist[, -the.edge.idx], node2)
  #nodes2 = setdiff(1:n, c(the.edge.idx, nodes1))

  # catf("Sel edge: %s", collapse(edgelist[, the.edge.idx]))
  # catf("Reachable from node %i: [%i] %s", node1, length(nodes1), collapse(nodes1))
  # catf("Reachable from node %i: [%i] %s", node2, length(nodes2), collapse(nodes2))
  # catf("Intersection: %i", length(intersect(nodes1, nodes2)))

  # now sample one node from each set and return build the corresponding edge

  # getRandomConnectingEdge = function(nodes1, nodes2) {
  #   adj.mat = grapherator::getAdjacencyMatrix(instance)

  #   # reduce to interesting submatrix
  #   adj.mat = adj.mat[nodes1, nodes2, drop = FALSE]

  #   # compute adjacency matrix
  #   adj.list = lapply(1:nrow(adj.mat), function(i) which(adj.mat[i, ]))

  #   # catf("Adj.mat: %s", collapse(dim(adj.mat)))
  #   # catf("nodes1: %i, nodes2: %i,  adj.list: %i", length(nodes1), length(nodes2), length(adj.list))

  #   #print(adj.mat)

  #   # which edges between nodes1 and nodes2 exist?
  #   adj.length = sapply(adj.list, length)
  #   idx.nonempty = which(adj.length != 0) # at least one must be

  #   # get random stuff
  #   idx.first = getRandom(idx.nonempty)
  #   idx.second = getRandom(adj.list[[idx.first]])

  #   # print(idx.nonempty)
  #   # print(adj.list)
  #   # print(idx.first)
  #   # print(idx.second)

  #   return(c(nodes1[idx.first], nodes2[idx.second]))
  # }

#   getRandomConnectingEdge = function(nodes1, nodes2) {
#     adj.mat = grapherator::getAdjacencyMatrix(instance)

#     # compute adjacency list
#     adj.list = lapply(1:nrow(adj.mat), function(i) which(adj.mat[i, ]))

#     # get adjacency list only for node set 1
#     adj.list2 = adj.list[nodes1]

#  #   catf("nodes1: %i, nodes2: %i,  adj.list: %i", length(nodes1), length(nodes2), length(adj.list))

#     # remove nodes within node set 1
#     adj.list2 = lapply(adj.list2, setdiff, nodes1)
# #    print(adj.list2)

#     adj.length = sapply(adj.list2, length)
#     idx.nonempty = which(adj.length > 0)

#     idx.first = getRandom(idx.nonempty)
#     idx.second = getRandom(adj.list2[[idx.first]])

#     # print(idx.nonempty)
#     # print(idx.first)
#     # print(idx.second)

#     return(c(nodes1[idx.first], idx.second))
#   }

  getRandomConnectingEdge = function(nodes1, nodes2) {
    adj.mat = grapherator::getAdjacencyMatrix(instance)

    # reduce to part of adjacency matrix
    adj.mat2 = adj.mat[nodes1, nodes2, drop = FALSE]
    #catf("Dim of subadjmat: %s", collapse(dim(adj.mat2)))

    #print(adj.mat2)

    # get matrix rows and cols
    pos.edges = which(adj.mat2, arr.ind = TRUE)
    #print(pos.edges)
    if (!is.matrix(pos.edges))
      pos.edges = matrix(pos.edges, nrow = 1L)

    #print(nrow(pos.edges))
    # sample an edge at random
    idx.newnodes = pos.edges[getRandom(1:nrow(pos.edges)), ]
    #catf("New edge idx: %s", collapse(idx.newnodes))

    return(c(nodes1[idx.newnodes[1L]], nodes2[idx.newnodes[2L]]))
  }

  # replace edge with new edge
  edgelist[, the.edge.idx] = getRandomConnectingEdge(nodes1, nodes2)
  #checkValidity(edgelist, instance)

  return(edgelist)
}

edgeExchange = function(edgelist, p = 1 / ncol(edgelist), instance = NULL) {
  # get number of tree edges
  m = ncol(edgelist)

  # replace each edge with prob p
  for (i in 1:m) {
    if (runif(1L) < p) {
      # sample random edge in tree
      edgelist = oneEdgeExchange(edgelist, i, instance = instance)
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


subgraphMST = function(edgelist, sigma, scalarize, instance) {
  m = ncol(edgelist)
  n.objectives = instance$n.weights
  nsel = sample(3:sigma, 1L)
  #catf("Selecting connected subgraph with >= %i nodes.", nsel)
  inds = 1:m
  # select random edge in tree as the starting point
  start = sample(inds, 1L)
  #catf("First edge: %i, (%i, %i)", start, edgelist[1, start], edgelist[2, start])
  sel.edges = start
  # the incident nodes of the edge determine the first
  sel.nodes = edgelist[, start]
  # walk through tree randomly
  #cur.node = sel.nodes[1L]
  # loop until we got a sufficiently big subtree
  while (length(sel.nodes) < nsel) {
    # check which edges are incident to selected nodes
    rr = apply(edgelist, 2L, function(x) any(x %in% sel.nodes))
    #catf("Edges incident to nodes %s: %s", collapse(sel.nodes), collapse(which(rr)))
    # filter already selected
    rr = setdiff(which(rr), sel.edges)
    #catf("Edges incident to nodes %s: %s", collapse(sel.nodes), collapse(rr))
    #FIXME: this can be made better
    sel.edges = sort(c(sel.edges, rr))
    sel.nodes = unique(as.integer(edgelist[, sel.edges]))
  }
  # now extract subgraph and apply Prim
  sel.nodes = sort(sel.nodes)
  #catf("Finally extracted %i nodes %s", length(sel.nodes), collapse(sel.nodes))
  dd = NULL
  if (!scalarize) {
    obj.idx = sample(1:n.objectives, 1L)
    dd = instance$weights[[obj.idx]][sel.nodes, sel.nodes]
  } else {
    lambdas = sampleWeights(n.objectives)
    dd = scalarizeWeights(instance$weights, lambdas)[sel.nodes, sel.nodes]
    # lambda = runif(1L)
    # dd = lambda * instance$weights[[1L]][sel.nodes, sel.nodes] + (1 - lambda) * instance$weights[[2L]][sel.nodes, sel.nodes]
  }

  #catf("submatrix:")
  #print(dd)
  # get result of PRIM
  mstres = vegan::spantree(d = dd)

  #print(mstres)
  #catf("N(sel.edges): %i", length(sel.edges))
  #catf("N(mstnodes):  %i", length(mstres$kid))
  repl.edges = matrix(
    #FIXME: is this correct?!?
    c(sel.nodes[2:(length(sel.edges) + 1L)],
      sel.nodes[mstres$kid]), byrow = TRUE, nrow = 2L)

    #catf("Replacing edges %s:", collapse(sel.edges))
    #print(edgelist[, sel.edges])
    #catf("... with ")
    #print(repl.edges)
  edgelist[, sel.edges] = repl.edges
  return(edgelist)
}

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
  mutator = function(ind, sigma = floor(ncol(ind) / 2), scalarize = FALSE, instance = NULL) {
    subgraphMST(ind, sigma, scalarize, instance)
  },
  supported = "custom"
)
