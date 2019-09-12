\name{DMR.plot}
\alias{DMR.plot}
\alias{plot}
\title{
  Plotting DMRs
}
\description{
  Plots an individual DMR (in context of possibly other DMRs) as found by \code{dmrcate}. 
  Heatmaps are shown as well as proximal coding regions, smoothed group means and chromosome ideogram.
}
\usage{
DMR.plot(ranges, dmr, CpGs, what=c("Beta", "M"), 
         arraytype=c("EPIC", "450K"), phen.col, 
         genome = c("hg19", "hg38", "mm10"), ...)
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
  \item{\dots}{
    Extra arguments passed to \code{Gviz:::plotTracks()}.
  }
}
\value{
  A plot to the current device.
}
\author{
  Aaron Statham <a.statham@garvan.org.au>, Tim J. Peters <t.peters@garvan.org.au>
}
\examples{

data(dmrcatedata)
myMs <- logit2(myBetas)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
patient <- factor(sub("-.*", "", colnames(myMs)))
type <- factor(sub(".*-", "", colnames(myMs)))
design <- model.matrix(~patient + type) 
myannotation <- cpg.annotate("array", myMs.noSNPs, what="M", arraytype = "450K",
                             analysis.type="differential", design=design, coef=39)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
groups <- c(Tumour="magenta", Normal="forestgreen")
cols <- groups[as.character(type)]
DMR.plot(ranges=results.ranges, dmr=1, CpGs=myBetas[,c(1:6, 38+(1:6))], what="Beta", arraytype = "450K", 
         phen.col=cols[c(1:6, 38+(1:6))], genome="hg19")

}
