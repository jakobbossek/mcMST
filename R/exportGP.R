exportGP = function(graph, file) {
  assertClass(graph, "mcGP")
  con = file(path.expand(file), "w")
  on.exit(close(con))

  n = graph$n.nodes
  m = if (!is.null(graph$adj.mat)) sum(graph$adj.mat) else n * n
  p = graph$n.weights

  meta = sprintf("%i,%i,%i,%i", n, m, graph$n.clusters, p)
  writeLines(meta, con)

  weight.types = collapse(graph$weight.types, ",")
  cat(weight.types, file = con, sep = "\n")

  node.types = collapse(graph$node.types, ",")
  cat(node.types, file = con, sep = "\n")

  write.table(graph$coordinates, file = con, append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)

  weights = matrix(NA, nrow = m, ncol = 2L + p)

  k = 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (!is.null(graph$adj.mat)) {
        if (graph$adj.mat[i, j] == 0L)
          next
      }
      weights[k, ] = c(i, j, sapply(graph$weights, function(weight.mat) weight.mat[i, j]))
      k = k + 1L
    }
  }
  write.table(weights, file = con, append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
  assert(k == m + 1)
  return(invisible(file))
}

importGP = function(file) {
  assertFile(file, access = "r")
  con = file(file, "r")
  on.exit(close(con))

  # generate bare graph
  g = mcGP(lower = 0, upper = 100)

  # import meta data
  meta = scan(file = con, what = integer(), n = 4L, sep = ",")
  g$n.nodes = n = meta[1L]
  g$n.clusters = meta[3L]
  g$n.weights = p = meta[4L]
  g$n.edges = m = meta[2L]

  # import weight and node types
  g$weight.types = scan(con, what = character(), n = p, sep = ",")
  g$node.types = strsplit(scan(con, what = character(), n = 1L), ",")[[1L]]

  # import coordinates
  g$coordinates = read.table(con, sep = ",", nrows = n, header = FALSE, stringsAsFactors = FALSE)

  # import edge->weight mapping
  ww = read.table(con, sep = ",", nrows = m, header = FALSE, stringsAsFactors = FALSE)

  # recreate weight matrices and adjacency matrix
  adj.mat = matrix(0, nrow = n, ncol = n)
  weights = vector(mode = "list", length = p)
  for (i in seq_len(p)) {
    weights[[i]] = matrix(Inf, nrow = n, ncol = n)
  }
  print(ww)
  for (k in seq_len(m)) {
    i = ww[k, 1L]
    j = ww[k, 2L]
    cur.weight = ww[k, 3:(2 + p)]
    for (l in seq_len(p)) {
      weights[[l]][i, j] = as.numeric(cur.weight[l])
    }
    adj.mat[i, j] = 1L
  }
  g$adj.mat = adj.mat
  g$weights = weights
  return(g)
}
