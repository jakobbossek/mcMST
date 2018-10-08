#' @export
readGP = function(x) {
  importer = new(GraphImporter)
  g = importer$importFromGrapheratorFile(x)
  g$saveVectorOfEdges()
  return(g)
}
