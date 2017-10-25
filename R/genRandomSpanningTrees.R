#' @title Generate a random spanning tree.
#'
#' @description Generate a random spanning tree of a graph
#' given the number of nodes of the problem instance.
#'
#' @template arg_n
#' @param type [\code{character(1)}]\cr
#'   String representing the desired format of the generated
#'   spanning tree. Possible values are \dQuote{pruefer} (Pruefer-code),
#'   \dQuote{edgelist} and \dQuote{charvec} (characteristic vector).
#'   Default is \dQuote{pruefer}.
#' @return [\code{integer} | \code{matrix(2, n)}] Return type depends on \code{type}.
#' @examples
#' genRandomSpanningTree(10)
#' genRandomSpanningTree(10, type = "edgelist")
#' @export
genRandomSpanningTree = function(n, type = "pruefer") {
  n = asInt(n, lower = 4L)
  assertChoice(type, choices = c("pruefer", "edgelist", "charvec"))
  pcode = sample(1:n, size = n - 2L, replace = TRUE)
  if (type == "pruefer")
    return(pcode)
  else if (type == "edgelist")
    return(prueferToEdgeList(pcode))
  else
    return(prueferToCharVec(pcode))
}

#' @title Generate a set of random spanning trees.
#'
#' @description Generate a set of random spanning trees of a graph
#' given the number of nodes of the problem instance.
#'
#' @inheritParams genRandomSpanningTree
#' @param m [\code{integer(1)}]\cr
#'   Number of random spanning trees to be generated.
#' @param simplify [\code{logical(1)}]\cr
#'   Should the result be simplified to a matrix if appropriate?
#'   Only relevant if \code{type} is either \dQuote{pruefer} or
#'   \dQuote{charvec}.
#'   Default is \code{TRUE}.
#' @return [\code{list} | \code{matrix}] Result type depends on \code{simplify}
#'   and \code{type}.
#' @examples
#' genRandomSpanningTrees(3, 10)
#' genRandomSpanningTrees(3, 10, simplify = FALSE)
#'
#' genRandomSpanningTrees(3, 10, type = "edgelist")
#' @export
genRandomSpanningTrees = function(m, n, type = "pruefer", simplify = TRUE) {
  m = asInt(m, lower = 1L)
  assertChoice(type, choices = c("pruefer", "edgelist", "charvec"))
  assertFlag(simplify)
  trees = lapply(seq_len(m), function(i) genRandomSpanningTree(n = n, type = type))
  if (simplify & (type %in% c("pruefer", "charvec")))
    return(do.call(rbind, trees))
  return(trees)
}
