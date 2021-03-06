#' Convert Pruefer code to edge list.
#'
#' @template arg_pcode
#' @return [\code{matrix(2, length(pcode) + 1)}] Edge list.
#' @examples
#' # here we generate a random Pruefer-code representing
#' # a random spanning tree of a graph with n = 10 nodes
#' pcode = sample(1:10, 8, replace = TRUE)
#' print(pcode)
#' edgelist = prueferToEdgeList(pcode)
#' print(edgelist)
#' @family transformation functions
#' @export
prueferToEdgeList = function(pcode) {
  n = length(pcode) + 2L
  # n = instance$n
  # if (length(pcode) != (n - 2))
  #   stopf("Prüfer-code of %i-vertex tree needs to have length %i, but has length %i!", n, n - 2, length(pcode))
  parents = integer(n - 1L)
  kids = integer(n - 1L)
  # now decode Prüfer-code

  # complement(P) is the set of all node numbers which are not in the code
  PC = setdiff(seq_len(n), pcode)
  # P is the original pcode
  P = pcode

  for (k in 1:(n - 2)) {
    j = which.min(PC)
    j1 = PC[j]
    kids[k] = j1
    parents[k] = P[k]
    PC = PC[-j]
    #catf("Removing %i from PC: ", j1)
    # add P[k] to PC if it does not exist in the remainder sequence of P
    #catf("checking if %i exists in P", P[k])
    if (k == (n - 2)) {
      #catf("NO")
      PC = c(PC, P[k])
      break
    }
    # if (is.na(!any(P[k] == P[(k+1):(n-2)]))) {
    #   print(k)
    #   print(P[k])
    #   print(pcode)
    #   print(P[(k+1):(n-2)])
    #   print(P)
    #   print(PC)
    # }
    if (!any(P[k] == P[(k+1):(n-2)])) {
      PC = c(PC, P[k])
    }
    #catf("PC: %s", collapse(PC))
    #catf("P: %s", collapse(P))
    #catf("-------")
  }
  #print(PC)

  # now add edge between last two remaining nodes in PC
  kids[n - 1L] = PC[1L]
  parents[n - 1L] = PC[2L]

  edges = rbind(kids, parents)

  return(edges)
}

#' Convert edge list to characteristic vector.
#'
#' @template arg_edgelist
#' @template arg_n
#' @template ret_charvec
#' @examples
#' # first we generate a small edge list by hand
#' # (assume the given graph has n = 4 nodes)
#' edgelist = matrix(c(1, 2, 2, 4, 3, 4), ncol = 3)
#' print(edgelist)
#' # next we transform the edge into
#' # a characteristic vector
#' cvec = edgeListToCharVec(edgelist, n = 4)
#' print(cvec)
#' @family transformation functions
#' @export
edgeListToCharVec = function(edgelist, n = NULL) {
  # number of nodes is |E| + 1
  n.edges = ncol(edgelist)
  if (is.null(n))
    n = n.edges + 1L

  mat = matrix(0, nrow = n, ncol = n)
  for (i in 1:n.edges) {
    tmp = sort(edgelist[, i])
    mat[tmp[1L], tmp[2L]] = 1L
  }
  # convert matrix colwise into a vector
  cv = as.integer(mat)
  stopifnot(sum(cv) == n.edges)
  return(cv)
}

#' Convert sequence of nodes to edge list.
#'
#' @param nodelist [\code{integer}]\cr
#'   Sequence of nodes.
#' @template ret_edgelist
#' @examples
#' # first generate a random permutation, e.g., representing
#' # a roundtrip tour in a graph
#' nodelist = sample(1:8)
#' # now convert into an edge list
#' nodelistToEdgelist(nodelist)
#' @family transformation functions
#' @export
nodelistToEdgelist = function(nodelist) {
  n = length(nodelist)
  edgelist = matrix(NA, nrow = 2L, ncol = n - 1L)
  for (i in 1:(n - 1L)) {
    edgelist[, i] = nodelist[i:(i + 1L)]
  }
  return(edgelist)
}

#' Convert characteristic vector to edge list.
#'
#' @template arg_charvec
#' @template ret_edgelist
#' @examples
#' # here we generate a random Pruefer-code representing
#' # a random spanning tree of a graph with n = 10 nodes
#' pcode = sample(1:10, 8, replace = TRUE)#'
#' edgelist = charVecToEdgelist(prueferToCharVec(pcode))
#' @family transformation functions
#' @export
charVecToEdgelist = function(charvec) {
  n = sqrt(length(charvec))
  mat = matrix(charvec, nrow = n, ncol = n, byrow = FALSE)
  edgelist = t(which(mat == 1L, arr.ind = TRUE))
  return(edgelist)
}

#' Convert Pruefer code to characteristic vector.
#'
#' @template arg_pcode
#' @template ret_charvec
#' @examples
#' # here we generate a random Pruefer-code representing
#' # a random spanning tree of a graph with n = 10 nodes
#' pcode = sample(1:10, 8, replace = TRUE)
#' print(pcode)
#' print(prueferToCharVec(pcode))
#' @family transformation functions
#' @export
prueferToCharVec = function(pcode) {
  edgeListToCharVec(prueferToEdgeList(pcode))
}

#' Convert permutation to edge list.
#'
#' @template arg_perm
#' @return [\code{matrix(2, length(perm))}] Edge list.
#' @examples
#' # first generate a random permutation, e.g., representing
#' # a roundtrip tour in a graph
#' perm = sample(1:10)
#' print(perm)
#' # now convert into an edge list
#' permutationToEdgelist(perm)
#' @family transformation functions
#' @export
permutationToEdgelist = function(perm) {
  n = length(perm)
  # close path
  perm = c(perm, perm[1L])
  edgelist = matrix(NA, nrow = 2L, ncol = n)
  for (i in seq_len(n)) {
    edgelist[, i] = perm[i:(i + 1L)]
  }
  return(edgelist)
}

#' Convert permutation to characteristic vector.
#'
#' @template arg_perm
#' @template arg_n
#' @template ret_charvec
#' @examples
#' # first generate a random permutation, e.g., representing
#' # a roundtrip tour in a graph
#' perm = sample(1:10)
#' print(perm)
#' # now convert into an edge list
#' permutationToCharVec(perm, n = 10)
#' @family transformation functions
#' @export
permutationToCharVec = function(perm, n) {
  edgeListToCharVec(permutationToEdgelist(perm), n)
}
