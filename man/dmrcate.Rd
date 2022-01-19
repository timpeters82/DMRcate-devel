\name{dmrcate}
\alias{dmrcate}
\title{
DMR identification
}
\description{
  The main function of this package. Computes a kernel estimate
  against a null comparison to identify significantly differentially (or
  variable) methylated regions.
}
\usage{
dmrcate(object, 
           lambda = 1000,
           C=NULL,
           pcutoff = "fdr", 
           consec = FALSE, 
           conseclambda = 10, 
           betacutoff = NULL,
           min.cpgs = 2
           ) 
}
\arguments{
  \item{object}{
    A \code{\link{CpGannotated-class}}, created from \code{\link{cpg.annotate}} or \code{\link{sequencing.annotate}}. 
  }
  \item{lambda}{
    Gaussian kernel bandwidth for smoothed-function estimation. Also informs DMR
    bookend definition; gaps >= \code{lambda} between significant CpG sites
    will be in separate DMRs. Support is truncated at 5*\code{lambda}. Default is 1000
    nucleotides. See details for further info.
  }
  \item{C}{
    Scaling factor for bandwidth. Gaussian kernel is calculated where 
    \code{lambda}/\code{C} = sigma. Empirical testing shows for both Illumina and bisulfite sequencing data that, when \code{lambda=1000}, near-optimal prediction of sequencing-derived DMRs is obtained when \code{C} is approximately 2, i.e. 1 standard deviation of Gaussian kernel = 500 base pairs. Cannot be < 0.2.
  }
  \item{pcutoff}{
    Threshold to determine DMRs. Default implies indexing at the rate of individually significant CpGs and can be set on the \code{\link{CpGannotated-class}} object using \code{\link{cpg.annotate}}, \code{\link{sequencing.annotate}} or \code{\link{changeFDR}}. \strong{Default highly recommended} unless you are comfortable with the risk of Type I error. If manually specified, this value will be set on the highly permissive kernel-smoothed FDR values.
  }
  \item{consec}{
    Use \code{DMRcate} in consecutive mode. Treats CpG sites as equally spaced.
  }
  \item{conseclambda}{
    Bandwidth in \emph{CpGs} (rather than nucleotides) to use when
    \code{consec=TRUE}. When specified the variable \code{lambda} simply
    becomes the minumum distance separating DMRs.
  }
  \item{betacutoff}{
    Optional filter; removes any region from the results where the absolute mean beta shift is less than the given value. Only available for Illumina array data and results produced from DSS::DMLtest().
  }
   \item{min.cpgs}{
    Minimum number of consecutive CpGs constituting a DMR.
  }
}
\details{

  The values of \code{lambda} and \code{C} should be chosen with care. For array data, we currently recommend that half a kilobase represent 1 standard deviation of support (\code{lambda=1000} and \code{C=2}). If \code{lambda} is too small or \code{C} too large then the kernel estimator will not have enough support to significantly differentiate the weighted estimate from the null distribution. If \code{lambda} is too large then \code{dmrcate} will report very long DMRs spanning multiple gene loci, and the large amount of support will likely give Type I errors. If you are concerned about Type I errors we \code{highly} recommend using the default value of \code{pcutoff}, although this will return no DMRs if no DM CpGs are returned by \code{limma/DSS} either.


}
\value{
  A DMResults object.
}

\references{

Peters, T. J., Buckley, M.J., Chen, Y., Smyth, G.K., Goodnow, C. C. and Clark, S. J. (2021). Calling differentially methylated regions from whole genome bisulphite sequencing with DMRcate. \emph{Nucleic Acids Research}, \bold{49}(19), e109.

Peters T.J., Buckley M.J., Statham, A., Pidsley R., Samaras K., Lord R.V., Clark S.J. and  Molloy P.L. \emph{De novo} identification of differentially methylated regions in the human genome. \emph{Epigenetics & Chromatin} 2015, \bold{8}:6,  doi:10.1186/1756-8935-8-6

Wand, M.P. & Jones, M.C. (1995) \emph{Kernel Smoothing}. Chapman & Hall.

Duong T. (2013) Local significant differences from nonparametric
two-sample tests. \emph{ Journal of Nonparametric Statistics}. 2013
\bold{25}(3), 635-645.

}
\author{
Tim J. Peters <t.peters@garvan.org.au>, Mike J. Buckley <Mike.Buckley@csiro.au>, Tim Triche Jr. <tim.triche@usc.edu>
}
\examples{
library(ExperimentHub)
library(limma)
eh <- ExperimentHub()
FlowSorted.Blood.EPIC <- eh[["EH1136"]]
tcell <- FlowSorted.Blood.EPIC[,colData(FlowSorted.Blood.EPIC)$CD4T==100 |
                                colData(FlowSorted.Blood.EPIC)$CD8T==100]
detP <- minfi::detectionP(tcell)
remove <- apply(detP, 1, function (x) any(x > 0.01))
tcell <- tcell[!remove,]
tcell <- minfi::preprocessFunnorm(tcell)
#Subset to chr2 only
tcell <- tcell[seqnames(tcell) == "chr2",]
tcellms <- minfi::getM(tcell)
tcellms.noSNPs <- rmSNPandCH(tcellms, dist=2, mafcut=0.05)
tcell$Replicate[tcell$Replicate==""] <- tcell$Sample_Name[tcell$Replicate==""]
tcellms.noSNPs <- avearrays(tcellms.noSNPs, tcell$Replicate)
tcell <- tcell[,!duplicated(tcell$Replicate)]
tcell <- tcell[rownames(tcellms.noSNPs),]
colnames(tcellms.noSNPs) <- colnames(tcell)
assays(tcell)[["M"]] <- tcellms.noSNPs
assays(tcell)[["Beta"]] <- minfi::ilogit2(tcellms.noSNPs)
type <- factor(tcell$CellType)
design <- model.matrix(~type) 
myannotation <- cpg.annotate("array", tcell, arraytype = "EPIC",
                             analysis.type="differential", design=design, coef=2)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
}


