#' @title Uniform mutation for Pruefer code representation.
#'
#' @description \code{mutUniform2} replaces each component of a Pruefer code of length n - 2
#' with probability \code{p} with a random node number between 1 and n.
#'
#' @param ind [\code{integer}]\cr
#'   Pruefer code.
#' @param p [\code{numeric(1)}]\cr
#'   Probability of mutation of each component of \code{ind}.
#'   Default is \code{1 / length(ind)}.
#' @return [\code{integer}] Mutated Pruefer code.
#' @export
mutUniform2 = makeMutator(
  mutator = function(ind, p = 1 / length(ind)) {
    n = length(ind)
    nrs = 1:(n + 2L)
    r = runif(n) < p
    if (sum(r) > 0)
      ind[r] = sample(nrs, sum(r), replace = TRUE)
    return(ind)
  },
  supported = "custom")

oneEdgeExchange = function(edgelist, edge.id) {
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

  getReachableNodes = function(edgelist, node) {
    nodes = node
    nodes2 = nodes
    while (length(nodes2) > 0L) {
      the.node = nodes2[1L]
      # check which edges are incident to current node
      is.adjacent = apply(edgelist, 2L, function(x) any(x == the.node))
      # get all unique nodes
      adjacent.nodes = unique(as.integer(edgelist[, is.adjacent]))
      nodes2 = setdiff(c(nodes2, adjacent.nodes), the.node)
      nodes = unique(c(nodes, adjacent.nodes))
      edgelist = edgelist[, !is.adjacent, drop = FALSE]
    }
    return(nodes)
  }

  # now for each end node get all reachable nodes via DFS
  nodes1 = getReachableNodes(edgelist[, -the.edge.idx], node1)
  nodes2 = getReachableNodes(edgelist[, -the.edge.idx], node2)

  #catf("Reachable from node %i: %s", node1, collapse(nodes1))
  #catf("Reachable from node %i: %s", node2, collapse(nodes2))

  # now sample one node from each set and return build the corresponding edge
  new.node1 = getRandom(nodes1)
  new.node2 = getRandom(nodes2)

  # replace edge with new edge
  edgelist[, the.edge.idx] = c(new.node1, new.node2)
  return(edgelist)
}

#' @title One-edge-exchange mutator for edge list representation of spanning trees.
#'
#' @description Each edge is replaced with another feasible edge with probability 1/m where
#' m is the number of edges, i.e., in expectation one edge is replaced by each
#' mutation operation.
#'
#' @param ind [\code{matrix(2, m)}]\cr
#'   Matrix of edges (each column is one edge).
#' @param p [\code{numeric(1)}]\cr
#'   Probability of edge exchange.
#'   Default is \code{1 / ncol(ind)}.
#' @return [\code{matrix(2, m)}] Mutated edge list.
#' @export
mutEdgeExchange = makeMutator(
  mutator = function(ind, p = 1 / ncol(ind)) {
    # get number of tree edges
    m = ncol(ind)

    # replace each edge with prob p
    for (i in seq_len(m)) {
      if (runif(1L) < p) {
        ind = oneEdgeExchange(ind, i)
      }
    }
    return(ind)
  },
  supported = "custom")

#' @title Subgraph-mutator for edge list representation.
#'
#' @description \code{mutSubgraphMST} selects a random edge e = (u, v) and traverses
#' the tree starting form u and v respectively until a connected subtree of at most
#' k edges is selected. Then the subtree is replaced with the optimal spanning tree
#' regarding objective 1 or 2 with equal probability.
#'
#' @param ind [\code{matrix(2, m)}]\cr
#'   Matrix of edges (each column is one edge).
#' @param sigma [\code{integer()}]\cr
#'   Upper bound for the size of the selected subtree.
#' @param instance [any]\cr
#'   Multi-objective graph problem instance.
#' @return [\code{matrix(2, m)}] Mutated edge list.
#' @export
mutSubgraphMST = makeMutator(
  mutator = function(ind, sigma = floor(ncol(ind) / 2), instance = NULL) {
    requirePackages("vegan", why = "rmoco::mutSubgraphMST")
    m = ncol(ind)
    nsel = sample(3:sigma, 1L)
    #catf("Selecting connected subgraph with >= %i nodes.", nsel)
    # select random edge in tree as the starting point
    start = sample(1:m, 1L)
    #catf("First edge: %i, (%i, %i)", start, ind[1, start], ind[2, start])
    sel.edges = start
    # the incident nodes of the edge determine the first
    sel.nodes = ind[, start]
    # walk through tree randomly
    #cur.node = sel.nodes[1L]
    inds = 1:m
    # loop until we got a sufficiently big subtree
    while (length(sel.nodes) < nsel) {
      # check which edges are incident to selected nodes
      rr = apply(ind, 2L, function(x) any(x %in% sel.nodes))
      #catf("Edges incident to nodes %s: %s", collapse(sel.nodes), collapse(which(rr)))
      # filter already selected
      rr = setdiff(which(rr), sel.edges)
      #catf("Edges incident to nodes %s: %s", collapse(sel.nodes), collapse(rr))
      #FIXME: this can be made better
      sel.edges = sort(c(sel.edges, rr))
      sel.nodes = unique(as.integer(ind[, sel.edges]))
    }
    # now extract subgraph and apply Prim
    sel.nodes = sort(sel.nodes)
    #catf("Finally extracted %i nodes %s", length(sel.nodes), collapse(sel.nodes))
    obj = if (runif(1L) < 0.5) 1L else 2L
    dd = if (obj == 1L) {
      #catf("first obj.")
      instance$w1[sel.nodes, sel.nodes]
    } else {
      instance$w2[sel.nodes, sel.nodes]
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
    #print(ind[, sel.edges])
    #catf("... with ")
    #print(repl.edges)
    ind[, sel.edges] = repl.edges
    return(ind)
  },
  supported = "custom"
)
