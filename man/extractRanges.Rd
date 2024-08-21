\name{extractRanges}
\alias{extractRanges}
\title{
Create a GRanges object from \code{\link{dmrcate}} output.
}
\description{
Takes a DMResults object and produces the corresponding GRanges object. 
}
\usage{
extractRanges(dmrcoutput, genome = c("hg19", "hg38", "mm10")) 
}

\arguments{
  \item{dmrcoutput}{
    A DMResults object.
}
  \item{genome}{
    Reference genome for annotating DMRs with promoter overlaps. 
    Can be one of \code{"hg19"}, \code{"hg38"} or \code{"mm10"}.
    Ranges are assumed to map to the reference stated; there is no liftover.
  }
}
\value{
A GRanges object.
}
\author{
Tim Triche Jr. <tim.triche@usc.edu>, Tim Peters <t.peters@garvan.org.au>
}
\examples{
library(AnnotationHub)
library(GenomicRanges)
ah <- AnnotationHub()
EPICv2manifest <- ah[["AH116484"]]
chr21probes <- rownames(EPICv2manifest)[EPICv2manifest$CHR=="chr21"]
coords <- EPICv2manifest[chr21probes, "MAPINFO"]
stats <- rt(length(chr21probes), 2)
pvals <- pt(-abs(stats), 100)
fdrs <- p.adjust(2*pvals, "BH")
annotated <- GRanges(rep("chr21", length(stats)), IRanges(coords, coords), stat = stats,
                     rawpval = pvals, diff = 0, ind.fdr = fdrs, is.sig = fdrs < 0.05)
names(annotated) <- chr21probes
myannotation <- new("CpGannotated", ranges=annotated)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
extractRanges(dmrcoutput, genome = "hg38")
}
