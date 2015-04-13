\name{DMRcate-package}
\alias{DMRcate-package}
\alias{DMRcate}
\docType{package}
\title{
Illumina 450K methylation array spatial analysis
}
\description{
  \emph{De novo} identification and extraction of differentially
  methylated regions (DMR) in the human genome using 450k array
  data. \code{DMRcate} extracts and annotates differentially methylated regions
  (DMRs) using an array-bias corrected smoothed estimate. Functions are
  provided for filtering probes possibly confounded by SNPs and
  cross-hybridisation. Includes bedGraph generation, GRanges generation and plotting functions.
}
\author{
Tim J. Peters <Tim.Peters@csiro.au>
}
\references{
Peters T.J., Buckley M.J., Statham, A., Pidsley R., Samaras K., Lord R.V., Clark S.J. and  Molloy P.L. \emph{De novo} identification of differentially methylated regions in the human genome. \emph{Epigenetics & Chromatin} 2015, \bold{8}:6,  doi:10.1186/1756-8935-8-6
}

\examples{
data(dmrcatedata)
myMs <- logit2(myBetas)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
patient <- factor(sub("-.*", "", colnames(myMs)))
type <- factor(sub(".*-", "", colnames(myMs)))
design <- model.matrix(~patient + type) 
myannotation <- cpg.annotate(myMs.noSNPs, analysis.type="differential",
    design=design, coef=39)
dmrcoutput <- dmrcate(myannotation, lambda=1000)
makeBedgraphs(dmrcoutput=dmrcoutput, betas=myBetas, samps=c(1,39))
phen.col <- c(rep("orange", 38), rep("blue", 38))
DMR.plot(dmrcoutput=dmrcoutput, dmr=1, betas=myBetas, phen.col=phen.col,
    pch=16, toscale=TRUE, plotmedians=TRUE)
}
