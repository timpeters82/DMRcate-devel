\name{DMRcate-package}
\alias{DMRcate-package}
\alias{DMRcate}
\docType{package}
\title{
DMR calling from bisulfite sequencing and Illumina array data
}
\description{
  \emph{De novo} identification and extraction of differentially
  methylated regions (DMRs) in the human genome using Illumin array and bisulfite sequencing
  data. \code{DMRcate} extracts and annotates differentially methylated regions
  (DMRs) using a kernel-smoothed estimate. Functions are
  provided for filtering probes possibly confounded by SNPs and
  cross-hybridisation. Includes GRanges generation and plotting functions.
}
\author{
Tim J. Peters <t.peters@garvan.org.au>
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
myannotation <- cpg.annotate("array", myMs.noSNPs, what="M", arraytype = "450K",
                             analysis.type="differential", design=design, coef=39)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
groups <- c(Tumour="magenta", Normal="forestgreen")
cols <- groups[as.character(type)]
DMR.plot(ranges=results.ranges, dmr=1, CpGs=myBetas, phen.col=cols, genome="hg19")
}
