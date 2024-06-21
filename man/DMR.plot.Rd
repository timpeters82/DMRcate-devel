\name{DMR.plot}
\alias{DMR.plot}
\alias{plot}
\title{
  Plotting DMRs
}
\description{
  Plots an individual DMR (in context of possibly other DMRs) as found by \code{dmrcate}. 
  Heatmaps are shown as well as proximal coding regions, smoothed methylation values 
  (with an option for smoothed group means) and chromosome ideogram.
}
\usage{
DMR.plot(ranges, 
         dmr, 
         CpGs, 
         what = c("Beta", "M"), 
         arraytype = c("EPICv2", "EPICv1", "450K"),
         phen.col,
         genome = c("hg19", "hg38", "mm10"),
         labels = names(ranges),
         flank = 5000,
         heatmap = TRUE, 
         extra.ranges = NULL, 
         extra.title = names(extra.ranges))
}
\arguments{
  \item{ranges}{
    A GRanges object (ostensibly created by \code{extractRanges()}) 
    describing DMR coordinates.
  }
  \item{dmr}{
    Index of \code{ranges} (one integer only) indicating which DMR to be
    plotted.
  }
  \item{CpGs}{
    Either:
    
    - A matrix of beta values for plotting, with unique Illumina probe IDs
    as rownames.
    
    - A GenomicRatioSet, annotated with the appropriate array and data types
    
    - A BSseq object containing per-CpG methylation and coverage counts for
    the samples to be plotted
  }
  \item{what}{
    Does \code{CpGs} (if a matrix) contain Beta or M-values? Not needed
    if object is a GenomicRatioSet or BSseq object.
  }
  \item{arraytype}{
    Is \code{CpGs} (if a matrix) sourced from EPIC or 450K data? Not needed
    if object is a GenomicRatioSet or BSseq object.
  }
  \item{phen.col}{
    Vector of colors denoting phenotypes of \emph{all} samples described in
    \code{CpGs}. See vignette for worked example.
  }
  \item{genome}{
    Reference genome for annotating DMRs. Can be one of \code{"hg19"},
    \code{"hg38"} or \code{"mm10"}
  }
  \item{labels}{
    Vector of DMR names to be displayed. Defaults to \code{names(ranges)}.
  }
  \item{flank}{
    Size, in base pairs, of the plotted region either side of the DMR. Cannot be less than 10bp or greater than 10kb.
  }
  \item{heatmap}{
    Should the heatmap be plotted? Default is \code{TRUE}, but \code{FALSE} is useful when plotting large numbers of samples.
  }
  \item{extra.ranges}{
    Optional GRanges object. Will plot any range overlapping a DMR..
  }
  \item{extra.title}{
    Vector of names for ranges from \code{extra.ranges}. Defaults to \code{names(extra.ranges)}.
  }
}
\value{
  A plot to the current device.
}
\author{
  Tim J. Peters <t.peters@garvan.org.au>, Aaron Statham <a.statham@garvan.org.au>, Braydon Meyer <b.meyer@garvan.org.au>
}
\examples{
library(GenomicRanges)
library(AnnotationHub)
ah <- AnnotationHub()
EPICv2manifest <- ah[["AH116484"]]
dmrranges <- GRanges("chr2:86787856-86793994")
probes <- EPICv2manifest$IlmnID[EPICv2manifest$CHR=="chr2" & 
                                EPICv2manifest$MAPINFO > 86770000 &
                                EPICv2manifest$MAPINFO < 86810000]
probes <- probes[order(EPICv2manifest[probes, "MAPINFO"])]
object <- minfi::logit2(matrix(rbeta(length(probes)*10, 3, 1), 
                               length(probes), 10))
rownames(object) <- probes
object[9:35, 6:10] <- minfi::logit2(matrix(rbeta(135, 1, 3), 27, 5))
cols <- c(rep("forestgreen", 5), rep("magenta", 5))
names(cols) <- rep(c("Ctrl", "Treat"), each=5)
DMR.plot(dmrranges, dmr = 1, CpGs=object, what = "M", arraytype="EPICv2", 
         phen.col=cols, genome="hg38")
}
