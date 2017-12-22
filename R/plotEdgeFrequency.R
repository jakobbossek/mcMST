#' @title Visualization of edge frequency among solution set.
#'
#' @description Given a list of graphs and a list of solutions (encoded as edge
#' lists) for each graph the function generates each one plot. This is a
#' 2D-scatterplot of edge weights of the graph. Size and colour of each point
#' indicate the number of solutions the edge is part of.
#'
#' @param graphs [\code{list(grapherator)}]\cr
#'   List of \code{\link[grapherator]{grapherator}} graphs.
#' @param approx.sets [\code{list(list(matrix))}]\cr
#'   List of approximations sets.
#' @param facet.args [\code{list}]\cr
#'   Further arguments passed down to \code{\link[ggplot2]{facet_wrap}}.
#'   Only relevant if \code{length(graphs) > 1}.
#' @param names [\code{character}]\cr
#'   Optional names of the graph instances. Used for facetting. Default is
#'   \dQuote{Problem_i} with i ranging from 1 to \code{length(graphs)}.
#' @return [\code{\link[ggplot2]{ggplot}}]
#' @examples
#' g = genRandomMCGP(50L)
#' res = mcMSTEmoaBG(mu = 10L, max.iter = 50, instance = g, scalarize = TRUE)
#' \dontrun{
#' plotEdgeFrequency(list(g), list(res$pareto.set))
#' }
#' @family result visualization
#' @export
plotEdgeFrequency = function(graphs, approx.sets, facet.args = list(), names = NULL) { # nocov start
  assertList(graphs, types = "grapherator", min.len = 1L)
  n = length(graphs)
  if (length(approx.sets) != n)
    stopf("plotEdgeFrequency: number of approx.sets must be equal to number of graphs.")
  assertList(approx.sets, types = "list")

  if (is.null(names))
    names = paste0("Problem_", 1:n)
  assertCharacter(names, len = n, any.missing = FALSE, all.missing = FALSE)

  # generate data
  ggdf = lapply(1:n, function(i) {
    approx.set = approx.sets[[i]]
    graph = graphs[[i]]
    cvs = lapply(approx.set, edgeListToCharVec, n = getNumberOfNodes(graph))
    cvs = do.call(rbind, cvs)

    # now generate scatter plots of edge weights indicating the frequency
    # of occurence in Pareto-optimal trees
    inAnyTree = colSums(cvs) > 0L
    # filter (i,i) edges
    idx.zero = which(!as.logical(getAdjacencyMatrix(graph)))
    #idx.zero = which(as.numeric(graph$weights[[1L]]) == 0)# | duplicated(as.numeric(inst$w1)))
    w1s = as.numeric(graph$weights[[1L]])[-idx.zero]
    w2s = as.numeric(graph$weights[[2L]])[-idx.zero]

    freq = colSums(cvs)[-idx.zero] / length(approx.set)

    data.frame(c1 = w1s, c2 = w2s, Frequency = freq, Problem = rep(names[i], length(freq)))
  })
  ggdf = do.call(rbind, ggdf)
  pl = ggplot2::ggplot(ggdf, ggplot2::aes_string(x = "c1", y = "c2", size = "Frequency", color = "Frequency"))
  pl = pl + ggplot2::geom_point()
  pl = pl + viridis::scale_colour_viridis()
  pl = pl + ggplot2::xlab(expression(c[1](e)))
  pl = pl + ggplot2::ylab(expression(c[2](e)))
  #title = sprintf("Instance: %s", inst$name)
  #pl = pl + ggtitle(title)
  pl = pl + ggplot2::guides(color = ggplot2::guide_legend(), size = ggplot2::guide_legend())
  pl = pl + ggplot2::theme(legend.position = "top")
  if (length(graphs) > 1L) {
    facet.args.defaults = list(facets = as.formula("~as.factor(Problem)"), nrow = 1L)
    facet.args = BBmisc::insert(facet.args.defaults, facet.args)
    pl = pl + do.call(ggplot2::facet_wrap, facet.args)
  }
  return(pl)
} # nocov end
