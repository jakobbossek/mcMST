#' Store an mcGP instance to file.
#'
#' @param instance [\code{mcGP}]\cr
#'   Problem instance.
#' @param file [\code{character(1)}]\cr
#'   Path to file.
#' @return [invisible(logical(1))]
saveMCGP = function(instance, file) {
  assertCharacter(file)

  n = instance$n
  #
  res = matrix("", nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      res[i, j] = paste0("(", instance$w1[i, j], ",", instance$w2[i, j], ")")
    }
  }

  #res = outer(1:n, 1:n, function(i, j) paste0("(", instance$w1[i, j], ", ", instance$w2[i, j], ")"))
  write.table(res, file = file, col.names = FALSE, row.names = FALSE, quote = FALSE)

  return(invisible(TRUE))
}
