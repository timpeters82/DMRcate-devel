#Author: tim.triche@usc.edu, April 29th, 2014

extractCoords <- function(xx) 
{ 
  if(length(xx) > 1) {
    t(sapply(xx, extractCoords))
  } else {
    coords <- strsplit(xx, ':')[[1]]
    chrom <- coords[1]
    chromStart <- strsplit(coords[2], '-')[[1]][1]
    chromEnd <- strsplit(coords[2], '-')[[1]][2]
    return(c(chrom=chrom, chromStart=chromStart, chromEnd=chromEnd))
  }
}
