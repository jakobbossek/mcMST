#' @export
readGP = function(x) {
  importer = new(GraphImporter)
  g = importer$importFromGrapheratorFile(x)
  g$saveVectorOfEdges()
  g$setEdgeProbabilities(rep(1, g$getE()))
  return(g)
}
