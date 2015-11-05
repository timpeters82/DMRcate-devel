extractRanges <-function(dmrcoutput, genome=c("hg19", "hg38", "mm10")) 
{
  stopifnot(is(dmrcoutput, "dmrcate.output"))
  coords <- extractCoords(dmrcoutput$results$coord)
  coords <- cbind(coords, dmrcoutput$results[, c("no.cpgs", "minfdr", "Stouffer", "maxbetafc", "meanbetafc")])
  coords$chromStart <- as.integer(as.character(coords$chromStart))
  coords$chromEnd <- as.integer(as.character(coords$chromEnd))
  ranges <- makeGRangesFromDataFrame(coords, keep.extra.columns = TRUE)
  switch(genome, 
    hg19={tx=tx.hg19},
    hg38={tx=tx.hg38},
    mm10={tx=tx.mm10}
  )
  promsidx <- as.data.frame(findOverlaps(ranges, promoters(tx, 2000, 2000)))
  proms <- tapply(promsidx$subjectHits, promsidx$queryHits, function(x) tx[x])
  ranges$overlapping.promoters <- NA
  promcount <- 1
  for(i in 1:length(ranges)){
    if(i %in% as.numeric(names(proms))){
      ranges$overlapping.promoters[i] <- paste(proms[[promcount]]$tx_name, collapse=", ")
      promcount <- promcount + 1
    }
  }
  ranges
}
