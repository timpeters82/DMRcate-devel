\name{dmrcate}
\alias{dmrcate}
\title{
DMR identification
}
\description{
  The main function of this package. Computes a kernel density estimate
  against a null comparison to identify significantly differentially (or
  variable, or hyper-) methylated regions in hg19.
}
\usage{
dmrcate(object, bw = 1000, p.adjust.method = "BH", pcutoff = 0.05,
    consec=FALSE, consecbw=10)
}
\arguments{
  \item{object}{
    A class of type "annot", created from \code{\link{annotate}}. 
  }
  \item{bw}{
    Gaussian kernel bandwidth for density estimation. Also informs DMR
    bookend definition; gaps >= \code{bw} between significant probes
    will be in separate DMRs since they do not lend support to each
    other. Support is truncated at 3.7 s.d. Default is 1000
    nucleotides. See details for further info.
  }
  \item{p.adjust.method}{
    Method for p-value adjustment from the significance test. Default is
    \code{"BH"} (Benjamini-Hochberg).
  }
  \item{pcutoff}{
    P-value cutoff to determine DMRs.
  }
  \item{consec}{
    Use DMRcate in consecutive probe mode. Treats CpG sites as equally spaced.
  }
  \item{consecbw}{
    Bandwidth in \emph{probes} (rather than nucleotides) to use when
    \code{consec=TRUE}. When specified the variable \code{bw} simply
    becomes the minumum distance separating DMRs.
  }
}
\details{

  The value of \code{bw} should be chosen with care. We recommend a
  bandwidth of 1 kilobase (\code{bw=1000}). If \code{bw} is too small
  then the kernel estimator will not have enough support to
  significantly differentiate the weighted estimate from the null
  distribution. If it is too large then \code{dmrcate} will report very
  long DMRs spanning multiple gene loci.

  Many gene loci have lengths reaching into the hundreds of thousands of
  base pairs, so it is quite possible that multiple signficant regions
  will have identical values in \code{results$gene_assoc}. This is fine;
  these regions are distinct in that they are at the very least
  \code{bw} nucleotides apart, and is preferable to attempting collapse
  into a super-DMR by increasing \code{bw}.

}
\value{
  A list of 2 data frames, \code{input} and \code{results}. \code{input}
  contains the contents of the "DMannot", "VMannot" or "HMannot" object,
  plus calculated p-values:
  \itemize{
    \item{ID}{As per annotation object input}
    \item{weights}{As per annotation object input}
    \item{CHR}{As per annotation object input}
    \item{pos}{As per annotation object input}
    \item{gene}{As per annotation object input}
    \item{group}{As per annotation object input}
    \item{raw}{Raw p-values from the significance test}
    \item{fdr}{Adjusted p-values from the significance test}
  }
  \code{results} contains an annotated data.frame of significant
  regions, ranked by \code{pval}:
  \itemize{
    \item{gene_assoc}{Complete list of gene loci overlapping the region,
      comma-separated}
    \item{group}{Complete list of gene annotations (e.g. TSS1500, 5'UTR
      etc.) overlapping the region, comma-separated}
    \item{hg19coords}{Coordinates of the significant region in
      hg19. IGV-friendly.}
    \item{no.probes}{Number of probes constituting the significant
      region. Tie-breaker when sorting probes by \code{pval}. A few
      regions may report \code{no.probes=1}, which may seem
      counter-intuitive, but this is only because the adjacent probes are
      just below the significance threshold. Unless \code{pcutoff} is
      highly conservative, it is unlikely that these regions will report
      at the head of the sorted list.}
    \item{pval}{Minimum adjusted p-value from the probes constituting the
      significant region.}
  }
}
\references{
Wand, M.P. & Jones, M.C. (1995) \emph{Kernel Smoothing}. Chapman & Hall.

Godtliebsen, F., Marron, J.S. & Chaudhuri, P. (2002) Significance in
scale space for bivariate density estimation. \emph{Journal of
Computational and Graphical Statistics}, \bold{11}, 1-22.

Duong T. (2013) Local significant differences from nonparametric
two-sample tests. \emph{ Journal of Nonparametric Statistics}. 2013
\bold{25}(3), 635-645.

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
myannotation <- annotate(myMs.noSNPs, analysis.type="differential",
    design=design, coef=39, diff.metric="FC", paired=TRUE, pcutoff=0.01)
dmrcoutput <- dmrcate(myannotation, bw=1000)
}
}


