\name{changeFDR}
\alias{changeFDR}
\title{
Change the individual CpG FDR thresholding for a CpGannotated object.
}
\description{
Takes a \code{\link{CpGannotated-class}} object and a specified FDR > 0 and < 1, and re-indexes the object in order to call DMRs at the specified rate.
}
\usage{
changeFDR(annot, FDR)
}
\arguments{
  \item{annot}{A \code{\link{CpGannotated-class}} object, created by \link{cpg.annotate} or \link{sequencing.annotate}.}
  \item{FDR}{The desired individual CpG FDR, which will index the rate at which DMRs are called.}
}
\details{
The number of CpG sites called as significant by this function will set the post-smoothing threshold for DMR constituents in \code{\link{dmrcate}}.
}
\value{
A re-indexed \code{\link{CpGannotated-class}} object.
}
\author{
Tim Peters <t.peters@garvan.org.au>
}
\examples{
library(GenomicRanges)
stats <- rt(1000, 2)
pvals <- pt(-abs(stats), 100)
fdrs <- p.adjust(2*pvals, "BH")
annotated <- GRanges(rep("chr1", 1000), IRanges(1:1000, 1:1000), stat = stats,
                     rawpval = pvals, diff = 0, ind.fdr = fdrs, is.sig = fdrs < 0.05)
names(annotated) <- paste0("CpG_", 1:1000)
myannotation <- new("CpGannotated", ranges=annotated)
changeFDR(myannotation, 0.1)
}

