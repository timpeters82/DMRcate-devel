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
library(ExperimentHub)
library(limma)
eh <- ExperimentHub()
FlowSorted.Blood.EPIC <- eh[["EH1136"]]
tcell <- FlowSorted.Blood.EPIC[,colData(FlowSorted.Blood.EPIC)$CD4T==100 |
                                colData(FlowSorted.Blood.EPIC)$CD8T==100]
detP <- detectionP(tcell)
remove <- apply(detP, 1, function (x) any(x > 0.01))
tcell <- tcell[!remove,]
tcell <- preprocessFunnorm(tcell)
#Subset to chr2 only
tcell <- tcell[seqnames(tcell) == "chr2",]
tcellms <- getM(tcell)
tcellms.noSNPs <- rmSNPandCH(tcellms, dist=2, mafcut=0.05)
tcell$Replicate[tcell$Replicate==""] <- tcell$Sample_Name[tcell$Replicate==""]
tcellms.noSNPs <- avearrays(tcellms.noSNPs, tcell$Replicate)
tcell <- tcell[,!duplicated(tcell$Replicate)]
tcell <- tcell[rownames(tcellms.noSNPs),]
assays(tcell, withDimnames=FALSE)[["M"]] <- tcellms.noSNPs
assays(tcell, withDimnames=FALSE)[["Beta"]] <- ilogit2(tcellms.noSNPs)
type <- factor(tcell$CellType)
design <- model.matrix(~type) 
myannotation <- cpg.annotate("array", tcell, arraytype = "EPIC",
                             analysis.type="differential", design=design, coef=2)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
groups <- c(CD8T="magenta", CD4T="forestgreen")
cols <- groups[as.character(type)]
DMR.plot(ranges=results.ranges, dmr=1, CpGs=getBeta(tcell), what="Beta", 
         arraytype = "EPIC", phen.col=cols, genome="hg19")
}
