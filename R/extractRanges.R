extractRanges <-function(dmrcoutput, genome=c("hg19", "hg38", "mm10"))
{
  env <- new.env(parent=emptyenv())
  data(dmrcatedata, envir=env)
  genome <- match.arg(genome)
  if(!is(dmrcoutput, "DMResults")){
    stop("Error: dmrcoutput is not a DMResults object. Please create one with dmrcate().")
  }
  coords <- extractCoords(dmrcoutput@coord)
  coords <- cbind(coords, dmrcoutput@no.cpgs, dmrcoutput@min_smoothed_fdr, dmrcoutput@Stouffer, 
                  dmrcoutput@HMFDR, dmrcoutput@Fisher, dmrcoutput@maxdiff, dmrcoutput@meandiff)
  coords$chromStart <- as.integer(as.character(coords$chromStart))
  coords$chromEnd <- as.integer(as.character(coords$chromEnd))
  ranges <- makeGRangesFromDataFrame(coords, keep.extra.columns = TRUE)
  switch(genome,
    hg19={grt=env$hg19.generanges},
    hg38={grt=env$hg38.generanges},
    mm10={grt=env$mm10.generanges}
  )
  genesidx <- as.data.frame(findOverlaps(ranges, grt))
  genesover <- tapply(genesidx$subjectHits, genesidx$queryHits, function(x) grt$symbol[x])
  op.A <- sapply(genesover, function(l) paste(l, collapse= ", "))
  name.A <- names(genesover)
  m.A <- as.numeric(name.A)
  M <- length(ranges)
  overlapping.genes <- rep(NA_character_, M)
  overlapping.genes[m.A] <- op.A
  ranges$overlapping.genes <- overlapping.genes
  colnames(values(ranges)) <- sub("dmrcoutput@", "", colnames(values(ranges)))
  ranges
}
