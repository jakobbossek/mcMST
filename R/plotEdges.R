#' @title Visualize edges common to several solutions.
#'
#' @description Given a list of characteristic vectors (graphs) the function plots
#' an embedding of the nodes in the Euclidean plane and depicts an edge if and
#' only if it is contained in at least one of the graphs. The edge thickness
#' indicates the number of graphs the edge is part of.
#'
#' @param x [\code{list}]\cr
#'   List of characteristic vectors.
#' @param n [\code{integer(1)}]\cr
#'   Number of nodes of the problem instance.
#'   Default is \code{sqrt(length(cv))} where cv is the first component of
#'   \code{x}.
#' @param normalize [\code{logical(1)}]\cr
#'   Shall edge frequencies be plotted? Default is code TRUE.
#' @param ... [any]\cr
#'   Further arguments passed down to \code{\link[qgraph]{qgraph}}.
#' @return Nothing
#' @family result visualization
#' @export
plotEdges = function(x, n = NULL, normalize = TRUE, ...) { # nocov start
  assertList(x, types = c("integer", "numeric"), min.len = 1L)
  if (is.null(n))
    n = sqrt(length(x))
  n = asInt(n, lower = 2)
  assertFlag(normalize)

  mat = matrix(0, nrow = n, ncol = n)
  for (i in 1:length(x)) {
    #print(x[i, ])
    mat[as.logical(x[[i]])] = mat[as.logical(x[[i]])] + 1L
  }

  if (normalize)
    mat = mat / length(x)

  qgraph.args = list(input = mat, layout = "spring", vsize = 1,
    directed = FALSE, edge.labels = TRUE,
    palette = "gray", theme = "gray", label.color = "white")
  qgraph.args = BBmisc::insert(qgraph.args, list(...))
  do.call(qgraph::qgraph, qgraph.args)
} # nocov end
