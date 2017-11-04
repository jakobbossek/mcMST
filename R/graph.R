#' @title Generate a bare multi-objective graph.
#'
#' @description This function generates a bare multi-objective graph. The generated
#' object does not contain nodes, edges or edge weights. It serves as a starting
#' point for a three step-approach of multi-objective graph problem construction:
#' 1) Add nodes respectively coordinates via \code{\link{addCoordinates}}, add edges
#' via \code{\link{addEdges}} and finally add edge weights with the function
#' \code{\link{addWeights}}.
#'
#' @param lower [\code{integer(1)}]\cr
#'   Lower bounds for coordinates.
#' @param upper [\code{integer(1)}]\cr
#'   Upper bounds for coordinates.
#' @template ret_mcGP
#' @family graph generators
#' @export
mcGP = function(lower, upper) {
  #n = asInt(n, lower = 2L)
  if (length(lower) == 1L)
    lower = rep(lower, 2L)
  if (length(upper) == 1L)
    upper = rep(upper, 2L)
  assertNumeric(lower, len = 2L, any.missing = FALSE, all.missing = FALSE)
  assertNumeric(upper, len = 2L, any.missing = FALSE, all.missing = FALSE)
  if (any(lower >= upper))
    stopf("mcGP: all elements of lower need to be stricly lower than the corresponding
      value in upper!")
  BBmisc::makeS3Obj(
    lower = lower,
    upper = upper,
    n.clusters = 0L,
    n.nodes = 0L,
    n.weights = 0L,
    weight.types = character(0L),
    node.types = character(0L),
    weights = list(),
    membership = NULL,
    coordinates = NULL,
    classes = "mcGP")
}

#' @export
print.mcGP = function(x, ...) {
  catf("MULTI-OBJECTIVE GRAPH PROBLEM")
  catf("Number of nodes: %i", x$n.nodes)
  if (x$n.clusters > 0L)
    catf("Number of clusters: %i", x$n.clusters)
  n.weights = length(x$weights)
  catf("Weights per edge: %i (%s)", n.weights, BBmisc::collapse(x$weight.types))
}

#' @title Coordinate generators.
#'
#' @description Functions for the placement of node coordinates in the
#' euclidean plane. Function \code{coordLHS} generates a space-filling
#' latin hypercube sample, \code{coordUniform} samples points from a
#' bivariate uniform distribution, \code{coordGrid} generates a regular
#' grid of points, \code{coordTriangular} generates a regular triangular
#' grid and \code{coordNormal} generates nodes on basis of a normal
#' distribution.
#'
#' @param n [\code{integer(1)}]\cr
#'   Number of points to generate.
#' @param lower [\code{numeric(2)}]\cr
#'   Minimal values for the first and second coordinates respectively.
#'   Default is 0.
#' @param upper [\code{numeric(2)}]\cr
#'   Maximal values for the first and second coordinates respectively.
#'   Default is 1.
#' @param method [\code{function}]\cr
#'   Function from package \pkg{lhs}.
#'   Default is \code{\link[lhs]{maximinLHS}}.
#' @param x.mean [\code{numeric}]\cr
#'   Mean value of normal distribution for x-value generation.
#'   Only relevant for \code{\link{coordNormal}}.
#' @param x.sd [\code{numeric}]\cr
#'   Standard deviation of normal distribution for x-value generation.
#'   Only relevant for \code{\link{coordNormal}}.
#' @param y.mean [\code{numeric}]\cr
#'   Mean value of normal distribution for y-value generation.
#'   Only relevant for \code{\link{coordNormal}}.
#' @param y.sd [\code{numeric}]\cr
#'   Standard deviation of normal distribution for y-value generation.
#'   Only relevant for \code{\link{coordNormal}}.
#' @return [\code{matrix(n, 2)}] Matrix of node coordinates.
#' @rdname coordGenerators
#' @name coordGenerators
#' @export
coordLHS = function(n, lower = 0, upper = 1, method = NULL) {
  if (is.null(method)) {
    requirePackages("lhs", why = "mcMST::coordLHS", default.method = "load")
    method = lhs::maximinLHS
  }

  coords = method(n, 2L)
  # stretch
  coords = lower + (upper - lower) * coords
  return(coords)
}

#' @export
#' @rdname coordGenerators
coordUniform = function(n, lower, upper) {
  coords = lapply(seq_len(2L), function(i) {
    runif(n, min = lower[i], max = upper[i])
  })
  coords = do.call(cbind, coords)
  return(coords)
}

#' @export
#' @rdname coordGenerators
coordTriangular = function(n, lower, upper) {
  m = sqrt(n)
  # determine offset of each second line
  d = (upper[1] - lower[1]) / (m - 1)
  d = d / 2
  print(d)
  print(m)
  offset = rep(c(rep(0, m), rep(d, m)), m)[1:n]
  print(offset)
  coords = coordGrid(n, lower, upper)
  coords[, 1L] = coords[, 1L] + offset
  return(coords)
}

#' @export
#' @rdname coordGenerators
coordGrid = function(n, lower, upper) {
  m = sqrt(n)
  x1 = seq(lower[1], upper[1], length.out = m)
  x2 = seq(lower[2], upper[2], length.out = m)
  coords = expand.grid(x1, x2)
  names(coords) = NULL
  coords = as.matrix(coords)
  return(coords)
}

#' @export
#' @rdname coordGenerators
coordNormal = function(n, lower, upper, x.mean, x.sd, y.mean, y.sd) {
  x1 = rnorm(n, x.mean, x.sd)
  x2 = rnorm(n, y.mean, y.sd)

  x1 = pmin(pmax(x1, lower[1L]), upper[1L])
  x2 = pmin(pmax(x2, lower[2L]), upper[2L])

  coords = cbind(x1, x2)
  return(coords)
}

#' @title Add node coordinates to graph.
#'
#' @description Places node coordinates in the two-dimensional euclidean plane.
#'
#' @template arg_mcGP
#' @param n [\code{integer}]\cr
#'   Number of coordinates to place. If \code{by.centers} is \code{FALSE} a single
#'   integer value is expected. Otherwise, a vector v may be passed. In this case
#'   v[i] coordinates are generated for each cluster. However, if a single value is
#'   passed and \code{by.center == TRUE}, each cluster is assigned the same number of
#'   coordinates.
#' @param coordinates [\code{matrix(n, 2)}]\cr
#'   Matrix of coordinates (each row is one point).
#'   Default is \code{NULL}. If this is set, setting of \code{generator}, \code{by.centers},
#'   and \code{par.fun} are ignored. This parameter is handy, if one wants to add
#'   coordinates by hand.
#'   Default is \code{NULL}.
#' @param generator [\code{function(n, ...)}]\cr
#'   Function used to generate coordinates. The generator needs to expect the number
#'   of points to generate as the first argument \code{n}. Additional control argument are
#'   possible.
#' @param by.centers [\code{logical(1)}]\cr
#'   Should coordinates be placed for each cluster center seperately? This enables
#'   geneation of clustered coordinates.
#'   Default is \code{FALSE}.
#' @param par.fun [\code{function(cc) | NULL}]\cr
#'   Optional function which is applied to each cluster center \code{cc} before the generation
#'   of coordinates in case \code{by.centers} is \code{TRUE}. This enables to specifically
#'   determine additional parameters for the \code{generator} for each cluster.
#' @param ... [any]\cr
#'   Furhter arguments passed down to \code{generator}.
#' @template ret_mcGP
#' @family graph generators
#' @export
addCoordinates = function(graph, n, generator, coordinates = NULL, by.centers = FALSE, par.fun = NULL, ...) {
  assertClass(graph, "mcGP")
  if (!is.null(coordinates)) {
    assertMatrix(coordinates, mode = "numeric", min.rows = 1L, ncols = 2L, any.missing = FALSE, all.missing = FALSE)
    n = nrow(coordinates)
  }
  n = asInt(n, lower = 1L)
  assertFlag(by.centers)
  assertFunction(par.fun, null.ok = TRUE)

  membership = NULL
  node.type = NULL

  # Helper function which aligns points with lower left point in [0,0].
  #
  # @param cluster.centers [matrix(2, n)]
  #   Matrix of city coordinates.
  # @return [matrix(2, n)]
  moveToOrigin = function(cluster.centers) {
    offset = abs(apply(cluster.centers, 2L, min))
    t(t(cluster.centers) + offset)
  }

  # if coordinates are passed, ignore the rest and add them
  if (!is.null(coordinates)) {
    coords = coordinates
    membership = if (graph$n.clusters == 0) rep(0, n) else rep(max(graph$membership) + 1L, n)
    node.type = "manual"
  }
  # if no two-phase approach simply delegate to coordinate generator
  else if (!by.centers) {
    coords = generator(n, lower = graph$lower, upper = graph$upper, ...)
    membership = if (graph$n.clusters == 0) rep(0, n) else rep(max(graph$membership) + 1L, n)
    node.type = "random"
  # otherwise use existing coordinates as cluster centers and place around them
  } else {
    # use current nodes as center coordinates
    center.coordinates = graph$coordinates
    if (graph$n.clusters > 0L)
      stopf("Currently one can add clusters only once!")
    nc = graph$n.nodes
    graph$n.clusters = nc
    graph$center.coordinates = center.coordinates
    # currently we allow only one "level" of clustering
    # thus we can set the membership here already
    graph$membership = 1:nc
    n = rep(n, nc)
    coords = lapply(seq_len(nc), function(i) {
      gen.args = list(n = n[i])
      # generate coordinates in origin
      if (!is.null(par.fun))
        gen.args = c(gen.args, par.fun(center.coordinates[i, ]))
      gen.args = c(gen.args, list(...))
      coords.cluster = do.call(generator, gen.args)

      #coords.cluster = moveToOrigin(coords.cluster)
      rects = apply(apply(coords.cluster, 2L, range), 2L, diff)
      cl.center = center.coordinates[i, ]
      # now move the way that centers are in fact centers
      #FIXME: ugly as hell
      coords.cluster = t(t(coords.cluster) + cl.center - rects / 2)
      return(coords.cluster)
    })
    # concatenate coordinates
    coords = do.call(rbind, coords)
    # assign membership (we know which cluster belongs to which center)
    membership = rep(1:nc, each = n[1L])
    node.type = "cluster"
  }
  # update meta data of graph
  graph$n.nodes = if (!is.null(graph$n.nodes)) graph$n.nodes + sum(n) else sum(n)
  graph$coordinates = if (!is.null(graph$coordinates)) rbind(graph$coordinates, coords) else coords
  graph$membership = if (!is.null(graph$membership)) c(graph$membership, membership)
  graph$node.types = c(graph$node.types, node.type)
  return(graph)
}

#' @title Define edges in multi-objective graph.
#'
#' @description By default \code{\link{addWeights}} generates n(n-1)/2 weights, i.e.,
#' the graph is assumed to be complete. This method allows to defne an adjacency
#' matrix to make the graph more sparse.
#'
#' @note Minimal implementation. No support so far.
#'
#' @template arg_mcGP
#' @param method [\code{function(...)}]\cr
#'   Method applied to \code{graph} in order to determine which edges to keep.
#'   Possible values are \dQuote{onion} or \dQuote{delauney}.
#' @family graph generators
#' @template ret_mcGP
addEdges = function(graph, method = "onion") { # nocov start
  assertClass(graph, "mcGP")
  assertChoice(method, choices = c("onion", "delauney"))

  adj.mat = matrix(0L, ncol = graph$n.nodes, nrow = graph$n.nodes)
  #diag(adj.mat) = 0

  # k-konvex hull approach
  #FIXME: links between node sets of onion layers
  if (method == "onion") {
    coordinates2 = graph$coordinates
    # which coordinates are already done?
    coords.done = rep(FALSE, graph$n.nodes)
    # indizes of coordinates not yet "onioned"
    # Needed as mapping for indices of coordinate matrix
    idx = which(!coords.done)
    n.edges = 0L
    while (TRUE) {
      # compute hull of remaining points
      ch = chull(coordinates2[!coords.done, , drop = FALSE])
      n.edges = n.edges + length(ch) - 1L
      # close hull
      ch = c(ch, ch[1L])
      # set edges
      for (i in (seq_along(ch) - 1L)) {
        adj.mat[idx[ch[i]], idx[ch[i + 1L]]] = 1L
        adj.mat[idx[ch[i + 1L]], idx[ch[i]]] = 1L
      }
      # update managing stuff
      coords.done[idx[ch]] = TRUE
      idx = which(!coords.done)
      if (sum(coords.done) == graph$n.nodes)
        break
    }
    print(n.edges)
  } else if (method == "delauney") {
    requirePackages("deldir", why = "mcMST:addEdges")
    # compute triangulation
    # The 5, 6 colums contains the indizes of the points
    dt = deldir(as.data.frame(graph$coordinates))$delsgs[, 5:6]
    for (i in seq_row(dt)) {
      adj.mat[dt[i, 1L], dt[i, 2L]] = 1L
      adj.mat[dt[i, 2L], dt[i, 1L]] = 1L
    }
  }
  graph$adj.mat = adj.mat
  return(graph)
} # nocov end

#' @title Add weights to a multi-objective graph.
#'
#' @description \code{addWeights} allows to generate edge weights for a multi-objective
#' graph instance. The weights can be generated on basis of the node coordinates (in this
#' case \code{\link[stats]{dist}} is applied with the cooresponding \code{method}).
#' Alternatively, all kinds of random weights can be generated.
#'
#' @template arg_mcGP
#' @param method [\code{character(1)}]\cr
#'   Method used to generate weights. Possible values are \dQuote{euclidean}, \dQuote{maximum},
#'   \dQuote{manhatten}, \dQuote{canberra}, \dQuote{binary}, \code{minkowski} or \code{random}.
#'   The latter generates (random) weights utilizing \code{weight.fun}. The remaining
#'   options are passed down to \code{\link[stats]{dist}}, i.e., weights are generated
#'   as distances between the node coordinates.
#' @param weights [\code{matrix}]\cr
#'   Square matrix of weights.
#'   If some weights are already assigned, pay attention to the correct dimensions.
#'   If this is passed all other arguments are ignored.
#'   Default is \code{NULL}.
#' @param weight.fun [\code{function(m, ...) | NULL}]\cr
#'   Function used to generate weights. The first arument needs to be number of weights
#'   to generate.
#' @param symmetric [\code{logical(1)}]\cr
#'   Should the weights be symmetric, i.e., w(i, j) = w(j, i) for each pair i, j of nodes?
#'   Default is \code{TRUE}.
#' @param to.int [\code{logical(1)}]\cr
#'   Should weights be rounded to integer?
#'   Default is \code{FALSE}.
#' @param rho [\code{numeric(1)}]\cr
#'   Correlation of edges weights for \code{method} \dQuote{correlated}.
#'   Default is \code{0.5}.
#' @param ... [any]\cr
#'   Additional arguments passed down to \code{weight.fun} or \code{\link[stats]{dist}}. See
#'   documentation of argument \code{method} for details.
#' @template ret_mcGP
#' @family graph generators
#' @export
addWeights = function(graph, method = "euclidean", weights = NULL, weight.fun = NULL, symmetric = TRUE, to.int = FALSE, rho = 0.5, ...) {
  assertClass(graph, "mcGP")
  assertChoice(method, choices = c("correlated", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "random"))
  assertFlag(to.int)
  assertNumber(rho, lower = -1, upper = 1)

  n.nodes = graph$n.nodes
  if (n.nodes == 0)
    stopf("addWeights: first place nodes/coordinates.")

  if (!is.null(weights))
    assertMatrix(weights, nrows = n.nodes, ncols = n.nodes, mode = "numeric")

  ws = graph$weights
  n.weights = if (is.null(ws)) 0L else length(ws)

  if (!is.null(weights)) {
    graph$weights[[n.weights + 1L]] = weights
    graph$n.weights = graph$n.weights + 1L
    graph$weight.types = c(graph$weight.types, "unknown")
  } else if (method == "correlated") {
    # get euclidean coordinates
    ww.euc = as.matrix(dist(graph$coordinates, method = "euclidean", ...))
    ww.euc.num = as.numeric(ww.euc)
    m = length(ww.euc.num)
    W = matrix(
      c(
        rep(1, m),
        ww.euc.num,
        runif(m, -1, 1)
      ),
    byrow = FALSE,
    ncol = 3L)

    # QR-decomposition
    Q = qr.Q(qr(W))

    T = matrix(c(1, rho, sqrt(1 - rho^2)), ncol = 3L)
    Y = T %*% t(Q)

    # normalize Y
    Y = (Y * graph$upper[1L])
    Y = Y + abs(min(Y)) + 10
    Y = matrix(Y, ncol = nrow(ww.euc))
    diag(Y) = 0

    if (!is.null(graph$adj.mat)) {
      ww.euc[graph$adj.mat == 0] = 10000 #FIXME: numeric infinity value
      Y[graph$adj.mat == 0] = 10000 #FIXME: numeric infinity value
    }

    graph$weights[[n.weights + 1L]] = ww.euc
    graph$weights[[n.weights + 2L]] = Y
    graph$n.weights = graph$n.weights + 2L
    graph$weight.types = c(graph$weight.types, c("euclidean", sprintf("%.2f-correlated", rho)))
  } else if (method != "random") {
    if (is.null(graph$coordinates))
      stopf("Method '%s' needs coordinates.", method)
    ww = as.matrix(dist(graph$coordinates, method = method, ...))

    if (!is.null(graph$adj.mat))
      ww[graph$adj.mat == 0] = 10000 #FIXME: numeric infinity value

    graph$weights[[n.weights + 1L]] = ww
    graph$n.weights = graph$n.weights + 1L
    graph$weight.types = c(graph$weight.types, "distance")
  } else {
    if (is.null(weight.fun))
      stopf("You need to pass a weight fun.")

    # always generate n^2 weights
    m = n.nodes * n.nodes

    ww = weight.fun(m, ...)

    #if (!is.null(adj.mat)) {
    ww = matrix(ww, ncol = n.nodes, nrow = n.nodes)
    diag(ww) = .0
    if (symmetric) {
      ww[lower.tri(ww)] = t(ww)[lower.tri(t(ww))]
    }

    if (!is.null(graph$adj.mat))
      ww[graph$adj.mat == 0] = 10000 #FIXME: numeric infinity value
    graph$weights[[n.weights + 1L]] = ww
    graph$n.weights = graph$n.weights + 1L
    graph$weight.types = c(graph$weight.types, "random")
  }

  return(graph)
}
