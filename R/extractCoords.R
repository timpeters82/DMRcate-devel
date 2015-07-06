extractCoords <- function(xx) 
{ 
  if(length(xx) > 1) {
    res <- data.frame(t(sapply(as.character(xx), extractCoords)))
    DataFrame(chrom=res$chrom, 
              chromStart=as.numeric(as.character(res$chromStart)), 
              chromEnd=as.numeric(as.character(res$chromEnd)))
  } else {
    coords <- strsplit(as.character(xx), ':')[[1]]
    chrom <- coords[1]
    chromStart <- as.numeric(strsplit(coords[2], '-')[[1]][1])
    chromEnd <- as.numeric(strsplit(coords[2], '-')[[1]][2])
    return(c(chrom=chrom, chromStart=chromStart, chromEnd=chromEnd))
  }
}
