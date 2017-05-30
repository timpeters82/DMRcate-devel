extractRanges <-function(dmrcoutput, genome=c("hg19", "hg38", "mm10"))
{
  env <- new.env(parent=emptyenv())
  data(dmrcatedata, envir=env)
  genome <- match.arg(genome)
  stopifnot(is(dmrcoutput, "dmrcate.output"))
  coords <- extractCoords(dmrcoutput$results$coord)
  coords <- cbind(coords, dmrcoutput$results[, c("no.cpgs", "minfdr", "Stouffer", "maxbetafc", "meanbetafc")])
  coords$chromStart <- as.integer(as.character(coords$chromStart))
  coords$chromEnd <- as.integer(as.character(coords$chromEnd))
  ranges <- makeGRangesFromDataFrame(coords, keep.extra.columns = TRUE)
  switch(genome,
    hg19={tx=env$tx.hg19},
    hg38={tx=env$tx.hg38},
    mm10={tx=env$tx.mm10}
  )
  promsidx <- as.data.frame(findOverlaps(ranges, promoters(tx, 2000, 2000)))
  proms <- tapply(promsidx$subjectHits, promsidx$queryHits, function(x) tx[x])
  op.A <- sapply(proms, function(l) paste(l$tx_name, collapse= ", "))
  name.A <- names(proms)
  m.A <- as.numeric(name.A)
  M <- length(ranges)
  overlapping.promoters <- rep(NA_character_, M)
  overlapping.promoters[m.A] <- op.A
  ranges$overlapping.promoters <- overlapping.promoters
  ranges
}
