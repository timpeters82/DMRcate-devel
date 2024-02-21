\name{DMRcate-package}
\alias{DMRcate-package}
\alias{DMRcate}
\docType{package}
\title{
DMR calling from bisulfite sequencing and Illumina array data
}
\description{
  \emph{De novo} identification and extraction of differentially
  methylated regions (DMRs) in the human genome using Illumina array and bisulfite sequencing
  data. \code{DMRcate} extracts and annotates differentially methylated regions
  (DMRs) using a kernel-smoothed estimate. Functions are
  provided for filtering probes possibly confounded by SNPs and
  cross-hybridisation. Includes GRanges generation and plotting functions.
}
\author{
Tim J. Peters <t.peters@garvan.org.au>
}
\references{
Peters TJ, Buckley MJ, Statham A, Pidsley R, Samaras K, Lord RV, Clark SJ and  Molloy PL. \emph{De novo} identification of differentially methylated regions in the human genome. \emph{Epigenetics & Chromatin} 2015, \bold{8}:6,  doi:10.1186/1756-8935-8-6

Peters TJ, Buckley MJ, Chen Y, Smyth GK, Goodnow CC and Clark SJ. Calling differentially methylated regions from whole genome bisulphite sequencing with DMRcate. \emph{Nucleic Acids Research} 2021, \bold{49}(19):e109. doi:10.1093/nar/gkab637.

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
myannotation <- cpg.annotate("array", tcell, 
                             arraytype = "EPICv1", analysis.type="differential",
                             design=design, coef=2)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
groups <- c(CD8T="magenta", CD4T="forestgreen")
cols <- groups[as.character(type)]
DMR.plot(ranges=results.ranges, dmr=1, 
         CpGs=minfi::getBeta(tcell), what="Beta", 
         arraytype = "EPICv1", phen.col=cols, genome="hg19")
}
