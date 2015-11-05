\name{DMR.plot}
\alias{DMR.plot}
\alias{plot}
\title{
  Plotting DMRs
}
\description{
  Plots an individual DMR (in context of possibly other DMRs) as found by \code{dmrcate}. 
  Heatmap and mean methylation plots are shown as well as genomic coordinates 
  and proximal coding regions.
}
\usage{
DMR.plot(ranges, dmr, CpGs, phen.col, genome = c("hg19", "hg38", "mm10"), 
                     array.annotation = c(array="IlluminaHumanMethylation450k", 
                                          annotation="ilmn12.hg19"),
                     samps = NULL, ...)
}
\arguments{
  \item{ranges}{
    A GRanges object (ostensibly created by \code{extractRanges()} 
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
    - A GRanges object describing individual CpGs to be plotted, containing
      methylated reads and total coverage for each sample. Please see the 
      worked example in the vignette for the correct structure of this object.
  }
  \item{phen.col}{
    Vector of colors denoting phenotypes of \emph{all} samples described in
    \code{CpGs}. See vignette for worked example.
  }
  \item{genome}{
    Reference genome for annotating DMRs. Can be one of \code{"hg19"},
    \code{"hg38"} or \code{"mm10"}
  }
  \item{array.annotation}{
    A vector describing the type of annotation from which plots are
    derived. Identical context to \code{minfi}, i.e. \code{annotation <-
    annotation(minfiobject)} where \code{minfiobject} is a
    \code{[Genomic](Methyl|Ratio)Set)}. Default (ilmn12.hg19) is
    recommended for array workflows.
  }
  
  \item{samps}{
    Vector of samples to be plotted, corresponding to indices of \code{phen.col}. 
    Default is all samples plotted.
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
\dontrun{
data(dmrcatedata)
myMs <- logit2(myBetas)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
patient <- factor(sub("-.*", "", colnames(myMs)))
type <- factor(sub(".*-", "", colnames(myMs)))
design <- model.matrix(~patient + type) 
myannotation <- cpg.annotate("array", myMs.noSNPs, analysis.type="differential",
    design=design, coef=39)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
groups <- c(Tumour="magenta", Normal="forestgreen")
cols <- groups[as.character(type)]
samps <- c(1:6, 38+(1:6))
DMR.plot(ranges=results.ranges, dmr=1, CpGs=myBetas, phen.col=cols, genome="hg19", samps=samps)
}
}
