extractRanges <- function(dmrcoutput) {
  stopifnot(is(dmrcoutput, 'dmrcate.output'))
  coords <- extractCoords(dmrcoutput$results$hg19coord)
  coords <- cbind(coords, dmrcoutput$results[,c("gene_assoc", "group", "no.probes",
                                                "minpval", "meanpval", "maxbetafc")])
  coords$chromStart <- as.integer(as.character(coords$chromStart))
  coords$chromEnd <- as.integer(as.character(coords$chromEnd))
  makeGRangesFromDataFrame(coords, keep.extra.columns=TRUE)
}
