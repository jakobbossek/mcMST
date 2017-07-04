# @title Pruefer code to edge list
#
# @description Function to transform a length (n-2) Prüfer code into an edge list of length
# (n-1)
#
# An edge list is stored as a (2 x (n-1)) matrix (each column represents one edge (i,j)).
# @param pcode [integer]
#   Spanning tree as a Pruefer code.
# @return [matrix] Edgelist.
prueferToEdgelist = function(pcode) {
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

# Get an edge-list and convert it into a characteristic vector.
#
# @param edgelist [(2 x n) matrix]
#   Matrix of edges (each column is one edge (i, j)).
# @return [integer] Characteristic vector cv with cv[i] = 1 if the i-th edge is in the tree.
edgeListToCharVec = function(edgelist, n = NULL) {
  # number of nodes is |E| + 1
  if (is.null(n))
    n = ncol(edgelist) + 1L

  mat = matrix(0, nrow = n, ncol = n)
  for (i in 1:(n - 1L)) {
    tmp = sort(edgelist[, i])
    mat[tmp[1L], tmp[2L]] = 1L
  }
  # convert matrix colwise into a vector
  cv = as.integer(mat)
  stopifnot(sum(mat) == n - 1L)
  return(as.integer(mat))
}

prueferToCharVec = function(pcode) {
  edgeListToCharVec(prueferToEdgelist(pcode))
}

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

permutationToCharVec = function(perm, n) {
  edgeListToCharVec(permutationToEdgelist(perm), n)
}
