# Load mcMST problem instance from a file.
#
# @param file [character(1)]
#   Path to file.
# @return [mcMST]
loadMCGP = function(file) {
  assertCharacter(file)

  res = read.csv(file, header = FALSE, sep = " ")
  n = nrow(res)
  w1 = matrix(NA, ncol = n, nrow = n)
  w2 = matrix(NA, ncol = n, nrow = n)

  for (i in 1:n) {
    for (j in 1:n) {
      tmp = gsub("\\(|\\)", "", res[i, j])
      tmp = as.numeric(strsplit(tmp, ",")[[1L]])
      w1[i, j] = tmp[1L]
      w2[i, j] = tmp[2L]
    }
  }
  BBmisc::makeS3Obj(classes = "mcGP", name = basename(file), n = n, w1 = w1, w2 = w2)
}
