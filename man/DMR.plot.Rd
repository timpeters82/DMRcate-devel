\name{DMR.plot}
\alias{DMR.plot}
\alias{plot}
\title{
  Plotting DMRs
}
\description{
  Plots an individual DMR as found by \code{dmrcate}.
}
\usage{
DMR.plot(dmrcoutput, dmr, betas, phen.col,
         annotation=c(array="IlluminaHumanMethylation450k",
                      annotation="ilmn12.hg19"),
         samps=NULL, toscale=FALSE, plotmedians=FALSE, ...)
}
\arguments{
  \item{dmrcoutput}{
    An object of class \code{dmrcate.output}.
  }
  \item{dmr}{
    Row index of \code{dmrcoutput$results}. Indicates which DMR to be
    plotted. Will only plot regions constituted of 2 or more CpGs.
  }
  \item{betas}{
    Matrix of beta values for plotting, with unique Illumina probe IDs
    as rownames.
  }
  \item{phen.col}{
    Vector of colors denoting phenotypes. Should be length
    \code{ncol(betas)[samps]}.
  }
  \item{annotation}{
    A vector describing the type of annotation from which plots are
    derived. Identical context to \code{minfi}, i.e. \code{annotation <-
    annotation(minfiobject)} where \code{minfiobject} is a
    \code{[Genomic](Methyl|Ratio)Set)}. Default (ilmn12.hg) is
    recommended.
  }
  \item{samps}{
    Vector of samples to be plotted, corresponding to columns of
    \code{betas}. Default is all samples plotted.
  }
  \item{toscale}{
    \code{TRUE} denotes CpGs plotted to scale along the x-axis according
    to their genomic coordinates, \code{FALSE} denotes evenly spaced
    plotting. Default is \code{FALSE}.
  }
  \item{plotmedians}{
    Logical denoting whether group medians will be plotted. Groups are derived from \code{phen.col}.
  }
  \item{\dots}{
    Extra arguments passed to \code{plot}
  }
}
\value{

  A plot to the current device. Square points along the top correspond
  to gene annotation; colours are as follows:

  \bold{TSS1500}: Light green

  \bold{TSS200}: Dark green

  \bold{Gene Body}: Red

  \bold{1st Exon}: Magenta

  \bold{5'UTR}: Dark Blue

  \bold{3'UTR}: Cyan
}
\author{
  Tim J. Peters <Tim.Peters@csiro.au>
}
\examples{
\dontrun{
data(dmrcatedata)
myMs <- logit2(myBetas)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
patient <- factor(sub("-.*", "", colnames(myMs)))
type <- factor(sub(".*-", "", colnames(myMs)))
design <- model.matrix(~patient + type) 
myannotation <- cpg.annotate(myMs.noSNPs, analysis.type="differential",
    design=design, coef=39)
dmrcoutput <- dmrcate(myannotation, lambda=1000)
phen.col <- c(rep("orange", 38), rep("blue", 38))
DMR.plot(dmrcoutput=dmrcoutput, dmr=1, betas=myBetas, phen.col=phen.col,
    pch=16, toscale=TRUE, plotmedians=TRUE)
}
}
