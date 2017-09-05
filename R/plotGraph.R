#' @title Visualize bi-objective graph.
#'
#' @description Only applicable for bi-objective problems of class \code{mcGP}.
#' \code{plot.mcGP} generates a scatterplot of edge weights. If the nodes do
#' have coordinates, additionally a scatterplot of the nodes in the euclidean
#' plane is generated.
#'
#' @param show.cluster.centers [\code{logical(1)}]\cr
#'   Display cluster centers?
#'   Default is \code{TRUE}. This option is ignored silently if the instance is not clustered.
#' @param ... [any]\cr
#'   Not used at the moment.
#' @param x [\code{mcGP}]\cr
#'   Multi-objective graph problem.
#' @param y Not used at the moment.
#' @return [\code{list}] A list of \code{\link[ggplot2]{ggplot}} objects with components
#' \code{pl.weights} (scatterplot of edge weights) and eventually \code{pl.coords} (scatterplot of
#' nodes). The latter is \code{NULL}, if \code{graph} has no associated coordinates.
#' @export
plot.mcGP = function(x, y=NULL, show.cluster.centers = TRUE, ...) {  
  
  assertFlag(show.cluster.centers)

  pl.coords = NULL
  n.nodes = x$n.nodes
  n.clusters = x$n.clusters
  n.weights = x$n.weights
  if (n.weights > 2L)
    stopf("plot.mcGP: More than 2 weights are currently not supported.")
  if (!is.null(x$coordinates)) {
    dd = as.data.frame(x$coordinates)
    names(dd) = c("x1", "x2")
    if (!is.null(x$membership))
      dd$Cluster = as.factor(x$membership)
    pl.coords = ggplot2::ggplot(dd, aes_string(x = "x1", y = "x2"))
    pl.coords = if (n.clusters > 0L) pl.coords + 
      ggplot2::geom_point(aes_string(colour = "Cluster")) else pl.coords + ggplot2::geom_point()
    if (n.clusters > 0L & show.cluster.centers) {
      ddc = as.data.frame(x$center.coordinates)
      names(ddc) = c("x1", "x2")
      pl.coords = pl.coords + 
        ggplot2::geom_point(data = ddc, colour = "black", size = 3)
      pl.coords = pl.coords +
        ggplot2::geom_point(data = ddc, colour = "white", size = 2.2)
    }
    pl.coords = pl.coords + 
      ggplot2::ggtitle("Node coordinates", subtitle = sprintf("#nodes: %i, #clusters: %i", n.nodes, n.clusters))
    pl.coords = pl.coords + 
      ggplot2::xlab(expression(x[1])) + 
      ggplot2::ylab(expression(x[2]))
  }
  weights1 = x$weights[[1L]]
  weights2 = x$weights[[2L]]

  weights1 = as.numeric(weights1[upper.tri(weights1)])
  weights2 = as.numeric(weights2[upper.tri(weights2)])

  dd = data.frame(w1 = weights1, w2 = weights2)
  pl.weights = ggplot2::ggplot(dd, aes_string(x = "w1", y = "w2")) + ggplot2::geom_point()
  pl.weights = pl.weights + 
    ggplot2::ggtitle("Edge weights", subtitle = collapse(paste(names(dd), x$weight.types, sep = " : "), sep = ","))
  pl.weights = pl.weights + xlab(expression(w[1])) + ylab(expression(w[2]))

  
  
  return(list(pl.coords = pl.coords, pl.weights = pl.weights))
}
