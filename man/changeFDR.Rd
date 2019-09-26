\name{changeFDR}
\alias{changeFDR}
\title{
Change the individual CpG FDR thresholding for a CpGannotated object.
}
\description{
Takes a \code{\link{CpGannotated-class}} object and a specified FDR > 0 and < 1, and re-indexes the object in order to call DMRs at the specified rate.
}
\usage{
changeFDR(annot, FDR)
}
\arguments{
  \item{annot}{A \code{\link{CpGannotated-class}} object, created by \link{cpg.annotate} or \link{sequencing.annotate}.}
  \item{FDR}{The desired individual CpG FDR, which will index the rate at which DMRs are called.}
}
\details{
The number of CpG sites called as significant by this function will set the post-smoothing threshold for DMR constituents in \code{\link{dmrcate}}.
}
\value{
A re-indexed \code{\link{CpGannotated-class}} object.
}
\author{
Tim Peters <t.peters@garvan.org.au>
}
\examples{
library(ExperimentHub)
library(limma)
library(edgeR)
eh <- ExperimentHub()
bis_1072 <- eh[["EH1072"]]
pData(bis_1072) <- data.frame(replicate=gsub(".*-", "", colnames(bis_1072)),
                              tissue=substr(colnames(bis_1072), 1, 
                                            nchar(colnames(bis_1072))-3), 
                              row.names=colnames(bis_1072))
colData(bis_1072)$tissue <- gsub("-", "_", colData(bis_1072)$tissue)
bis_1072 <- renameSeqlevels(bis_1072, mapSeqlevels(seqlevels(bis_1072), "UCSC"))
bis_1072 <- bis_1072[seqnames(bis_1072)=="chr2",]
tissue <- factor(pData(bis_1072)$tissue)
tissue <- relevel(tissue, "Liver_Treg")
design <- model.matrix(~tissue)
colnames(design) <- gsub("tissue|\\(|\\)", "", colnames(design))
rownames(design) <- colnames(bis_1072)
methdesign <- modelMatrixMeth(design)
cont.mat <- makeContrasts(treg_vs_tcon=Lymph_N_Treg-Lymph_N_Tcon,
                          fat_vs_ln=Fat_Treg-Lymph_N_Treg,
                          skin_vs_ln=Skin_Treg-Lymph_N_Treg,
                          fat_vs_skin=Fat_Treg-Skin_Treg,
                          levels=methdesign)
seq_annot <- sequencing.annotate(bis_1072, methdesign, all.cov = TRUE, 
                                   contrasts = TRUE, cont.matrix = cont.mat, 
                                   coef = "fat_vs_skin", fdr=0.05)
seq_annot <- changeFDR(seq_annot, 0.25)
}

