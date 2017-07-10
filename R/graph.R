#' @title Generate a bare multi-objective graph.
#'
#' @description This function generates a bare multi-objective weights. The generated
#' object does not contain nodes, edges or edge weights. It serves as a starting
#' point for the step-by-step construction of multi-objective graph problem.
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
    weights = list(),
    membership = NULL,
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
#' bivariate uniform distribution and \code{coordGrid} generates a regular
#' grid of points.
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
coordGrid = function(n, lower, upper) {
  m = sqrt(n)
  x1 = seq(lower[1], upper[2], length.out = m)
  x2 = seq(lower[2], upper[2], length.out = m)
  coords = expand.grid(x1, x2)
  names(coords) = NULL
  coords = as.matrix(coords)
  return(coords)
}

#' @title Add cluster centers to graph.
#'
#' @description Places \code{n.centers} cluster centers in the two-dimensional
#' euclidean plane by means of a \code{generator}, e.g., by Latin-Hypercube-Sampling (LHS).
#'
#' @template arg_mcGP
#' @param n.centers [\code{integer(1)}]\cr
#'   Number of cluster centers.
#' @param center.coordinates [\code{matrix(n, 2)}]\cr
#'   Matrix of center coordinates (each row is one point).
#'   Default is \code{NULL}. If this is set, \code{n.centers} and \code{generator}
#'   are ignored.
#' @param generator [\code{function(n, ...)}]\cr
#'   Function used to generate cluster centers. The generator needs to expect the number
#'   of points to generate as the first argument \code{n}. Additional control argument are
#'   possible.
#' @param ... [any]\cr
#'   Additional arguments passed down to \code{generator}.
#' @template ret_mcGP
#' @family graph generators
#' @export
addCenters = function(graph, n.centers = NULL, center.coordinates = NULL, generator = NULL, ...) {
  # sanity checks
  assertClass(graph, "mcGP")
  if (!is.null(n.centers))
    n.centers = asInt(n.centers, lower = 2L)
  assertFunction(generator, null.ok = TRUE)

  # more sanity checks
  if (!is.null(graph$coordinates))
    stopf("Graph already has coordinates! Place centers before coordinates.")
  if (!is.null(graph$center.coordinates))
    stopf("Cluster centers already placed.")
  if (is.null(n.centers) & is.null(center.coordinates))
    stopf("At least one of n.centers and center.coordinates must not be NULL.")

  if (is.null(center.coordinates))
    center.coordinates = generator(n.centers, lower = graph$lower, upper = graph$upper, ...)

  # generate cluster centers
  graph$center.coordinates = center.coordinates
  graph$n.clusters = nrow(center.coordinates)
  if (!("mcGP_clustered" %in% class(graph)))
    graph = addClasses(graph, "mcGP_clustered")
  return(graph)
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
addCoordinates = function(graph, n, generator, by.centers = FALSE, par.fun = NULL, ...) {
  assertClass(graph, "mcGP")
  n = asInteger(n, lower = 1L, any.missing = FALSE, all.missing = FALSE)
  if (length(n) > 1L & !by.centers)
    stopf("addCoordinates: vector of length > 1 for n only allowed if by.centers = TRUE.")
  assertFlag(by.centers)
  assertFunction(par.fun, null.ok = TRUE)

  membership = NULL

  # Helper function which aligns points with lower left point in [0,0].
  #
  # @param cluster.centers [matrix(2, n)]
  #   Matrix of city coordinates.
  # @return [matrix(2, n)]
  moveToOrigin = function(cluster.centers) {
    offset = abs(apply(cluster.centers, 2L, min))
    t(t(cluster.centers) + offset)
  }

  if (!by.centers) {
    coords = generator(n, lower = graph$lower, upper = graph$upper, ...)
  } else {
    nc = graph$n.clusters
    if (length(n) == nc) {
      n.per.cluster2 = n
    } else {
      n.per.cluster = floor(n / nc)
      n.per.cluster2 = rep(n.per.cluster, nc)
      if (nc * n.per.cluster != n) {
        idx = sample(seq_len(nc), 1L)
        n.per.cluster2[idx] = n.per.cluster2[idx] + 1L
      }
    }
    coords = lapply(seq_len(nc), function(i) {
      gen.args = list(n = n.per.cluster2[i])
      # generate coordinates in origin
      if (!is.null(par.fun))
        gen.args = c(gen.args, par.fun(graph$center.coordinates[i, ]))
      gen.args = c(gen.args, list(...))
      coords.cluster = do.call(generator, gen.args)

      coords.cluster = moveToOrigin(coords.cluster)
      rects = apply(apply(coords.cluster, 2L, range), 2L, diff)
      cl.center = graph$center.coordinates[i, ]
      # now move the way that centers are in fact centers
      #FIXME: ugly as hell
      coords.cluster = t(t(coords.cluster) + cl.center - rects / 2)
      return(coords.cluster)
    })
    # concatenate coordinates
    coords = do.call(rbind, coords)
    # assign membership (we know which cluster belongs to which center)
    membership = rep(1:nc, n.per.cluster2)
  }
  # update meta data of graph
  graph$n.nodes = if (!is.null(graph$n.nodes)) graph$n.nodes + sum(n) else sum(n)
  graph$coordinates = if (!is.null(graph$coordinates)) rbind(graph$coordinates, coords) else coords
  graph$membership = if (!is.null(graph$membership)) c(graph$membership, if (by.centers) membership else rep(0, nrow(coords))) else membership
  return(graph)
}

# @title Define edges in multi-objective graph.
#
# @description By default \code{\link{addWeights}} generates n(n-1)/2 weights, i.e.,
# the graph is assumed to be complete. This method allows to defne an adjacency
# matrix to make the graph more sparse.
#
# @note Minimal implementation. No support so far.
#
# @template arg_mcGP
# @param method [\code{function(...)}]\cr
#   Method applied to \code{graph} in order to determine which edges to keep.
# @family graph generators
# @template ret_mcGP
addEdges = function(graph, method = NULL) { # nocov start
  assertClass(graph, "mcGP")
  assertFunction(method, null.ok = TRUE)
  adj.mat = matrix(1L, ncol = graph$n.nodes, nrow = graph$n.nodes)
  diag(adj.mat) = 0
  graph$adj.mat = adj.mat
  return(graph)
} # nocov end

# edgesGrid = function(coordinates) {

# }

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
#' @param n [\code{integer(1)}]\cr
#'   Number of nodes. This is required only if there are no coordinates or no weights
#'   until now.
#'   Default is \code{NULL}, i.e., the number of nodes is extracted from \code{graph}.
#' @param symmetric [\code{logical(1)}]\cr
#'   Should the weights be symmetric, i.e., w(i, j) = w(j, i) for each pair i, j of nodes?
#'   Default is \code{TRUE}.
#' @param ... [any]\cr
#'   Additional arguments passed down to \code{weight.fun} or \code{\link[stats]{dist}}. See
#'   documentation of argument \code{method} for details.
#' @template ret_mcGP
#' @family graph generators
#' @export
addWeights = function(graph, method = "euclidean", weights = NULL, weight.fun = NULL, n = NULL, symmetric = TRUE, ...) {
  assertClass(graph, "mcGP")
  assertChoice(method, choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "random"))

  n.nodes = graph$n.nodes
  # check if graph has coordinates or weights already
  if (n.nodes == 0L) {
    if (!is.null(weights))
      n.nodes = ncol(weights)
    else
      n.nodes = n
    if (is.null(n.nodes))
      stopf("addWeights: number of nodes unknown. Please pass parameter 'n'.")
  }

  if (!is.null(weights))
    assertMatrix(weights, nrows = n.nodes, ncols = n.nodes, mode = "numeric")

  ws = graph$weights
  n.weights = if (is.null(ws)) 0L else length(ws)

  if (!is.null(weights)) {
    ww = weights
  } else if (method != "random") {
    if (is.null(graph$coordinates))
      stopf("Method '%s' needs coordinates.", method)
    ww = as.matrix(dist(graph$coordinates, method = method, ...))
    # if not all edges exist, set the remaining to infinifty
    # but keep zero distances on the diagonal
    if (!is.null(graph$adj.mat)) { # nocov start
      ww[graph$adj.mat == 0] = Inf
      diag(ww) = 0
    } # nocov end
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
      ww[graph$adj.mat == 0] = Inf #FIXME: numeric infinity value
  }
  graph$weights[[n.weights + 1L]] = ww
  wtype = if (!is.null(weights)) "manual" else ifelse(method != "random", "distance", "random")
  graph$weight.types = c(graph$weight.types, wtype)
  graph$n.weights = graph$n.weights + 1L

  if (graph$n.nodes == 0L)
    graph$n.nodes = graph$n.nodes + n.nodes
  return(graph)
}
