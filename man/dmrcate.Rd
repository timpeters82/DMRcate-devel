\name{dmrcate}
\alias{dmrcate}
\title{
DMR identification
}
\description{
  The main function of this package. Computes a kernel estimate
  against a null comparison to identify significantly differentially (or
  variable) methylated regions in hg19.
}
\usage{
dmrcate(object, 
        lambda = 1000,
        C=2,
        p.adjust.method = "BH", 
        pcutoff = "limma", 
        consec = FALSE, 
        conseclambda = 10, 
        betacutoff = NULL
        ) 
}
\arguments{
  \item{object}{
    A class of type "annot", created from \code{\link{cpg.annotate}}. 
  }
  \item{lambda}{
    Gaussian kernel bandwidth for smoothed-function estimation. Also informs DMR
    bookend definition; gaps >= \code{lambda} between significant probes
    will be in separate DMRs. Support is truncated at 5*\code{lambda}. Default is 1000
    nucleotides. See details for further info.
  }
  \item{C}{
    Scaling factor for bandwidth. Gaussian kernel is calculated where 
    \code{lambda}/\code{C} = sigma. Empirical testing shows that when \code{lambda=1000}, near-optimal prediction of sequencing-derived DMRs is obtained when \code{C} is approximately 2, i.e. 1 standard deviation of Gaussian kernel = 500 base pairs. Cannot be < 0.2.
  }
  \item{p.adjust.method}{
    Method for \emph{p}-value adjustment from the significance test. Default is
    \code{"BH"} (Benjamini-Hochberg).
  }
  \item{pcutoff}{
    \emph{p}-value cutoff to determine DMRs. Default is automatically determined by the number of significant probes returned by \code{limma} for that contrast, but can be set manually with a numeric value.
  }
  \item{consec}{
    Use \code{DMRcate} in consecutive probe mode. Treats CpG sites as equally spaced.
  }
  \item{conseclambda}{
    Bandwidth in \emph{probes} (rather than nucleotides) to use when
    \code{consec=TRUE}. When specified the variable \code{lambda} simply
    becomes the minumum distance separating DMRs.
  }
  \item{betacutoff}{
    Optional filter; removes any region from the results that does not have at least one         CpG site with a beta fold change exceeding this value.
  }
}
\details{

  The values of \code{lambda} and \code{C} should be chosen with care. We recommend that half a kilobase represent 1 standard deviation of support (\code{lambda=1000} and \code{C=2}). If \code{lambda} is too small or \code{C} too large then the kernel estimator will not have enough support to significantly differentiate the weighted estimate from the null distribution. If \code{lambda} is too large then \code{dmrcate} will report very long DMRs spanning multiple gene loci, and the large amount of support will likely give Type I errors. If you are concerned about Type I errors we recommend using the default value of \code{pcutoff}, although this will return no DMRs if no DM probes are returned by \code{limma} either.

  Many gene loci have lengths reaching into the hundreds of thousands of base pairs, so it is quite possible that multiple signficant regions will have identical values in \code{results$gene_assoc}. This is fine; these regions are distinct in that they are at the very least \code{lambda} nucleotides apart, and is preferable to attempting collapse
  into a super-DMR by increasing \code{lambda}.

}
\value{
  A list containing 2 data frames (\code{input} and \code{results}) and a numeric value (\code{cutoff}). \code{input} contains the contents of the \code{annot} object, plus calculated \emph{p}-values:
  \itemize{
    \item \code{ID}: As per annotation object input
    \item \code{weights}: As per annotation object input
    \item \code{CHR}: As per annotation object input
    \item \code{pos}: As per annotation object input
    \item \code{gene}: As per annotation object input
    \item \code{group}: As per annotation object input
    \item \code{betafc}: As per annotation object input
    \item \code{raw}: Raw \emph{p}-values from the significance test
    \item \code{fdr}: Adjusted \emph{p}-values from the significance test
  }
  \code{results} contains an annotated data.frame of significant
  regions, ranked by \code{minpval}:
  \itemize{
    \item \code{gene_assoc}: Complete list of gene loci overlapping the region,
      comma-separated
    \item \code{group}: Complete list of gene annotations (e.g. TSS1500, 5'UTR
      etc.) overlapping the region, comma-separated
    \item \code{hg19coords}: Coordinates of the significant region in
      hg19. IGV-friendly.
    \item \code{no.probes}: Number of probes constituting the significant
      region. Tie-breaker when sorting probes by \code{minpval}. A few
      regions may report \code{no.probes=1}, which may seem
      counter-intuitive, but this is only because the adjacent probes are either
      just below the significance threshold, or it is a highly DM probe in a sparse region. Unless \code{pcutoff} is highly conservative, it is unlikely that these regions will report at the head of the sorted list.
    \item \code{minpval}: Minimum adjusted \emph{p}-value from the probes constituting the
      significant region.
    \item \code{meanpval}: Mean adjusted \emph{p}-value from the probes constituting the
      significant region.
    \item \code{maxbetafc}: Maximum absolute beta fold change within the region
  }
  \code{cutoff} is the signficance \emph{p}-value cutoff provided in the call to \code{dmrcate}.
}

\references{

Peters T.J., Buckley M.J., Statham, A., Pidsley R., Samaras K., Lord R.V., Clark S.J. and  Molloy P.L. \emph{De novo} identification of differentially methylated regions in the human genome. \emph{Epigenetics & Chromatin} 2015, \bold{8}:6,  doi:10.1186/1756-8935-8-6

Wand, M.P. & Jones, M.C. (1995) \emph{Kernel Smoothing}. Chapman & Hall.

Duong T. (2013) Local significant differences from nonparametric
two-sample tests. \emph{ Journal of Nonparametric Statistics}. 2013
\bold{25}(3), 635-645.

}
\author{
Tim J. Peters <Tim.Peters@csiro.au>, Mike J. Buckley <Mike.Buckley@csiro.au>, Tim Triche Jr. <tim.triche@usc.edu>
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
}
}


