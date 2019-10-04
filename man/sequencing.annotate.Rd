\name{sequencing.annotate}
\alias{sequencing.annotate}
\title{
Annotate a bisulfite sequencing experiment (WGBS or RRBS) with probe weights and chromosomal position.
}
\description{
  Either:
  - Annotate a BSseq object with chromosome position and test statistic, or
  - Parse output from \code{DSS::DMLtest()} or \code{DSS::DMLtest.multiFactor()} into a CpGannotated object.
}
\usage{

sequencing.annotate(obj, methdesign, all.cov=FALSE, contrasts = FALSE, 
                                cont.matrix = NULL, fdr = 0.05, coef, ...) 
}
\arguments{
  \item{obj}{
    A BSseq object or data.frame output from \code{DSS::DMLtest()} or \code{DSS::DMLtest.multiFactor()}.
  }
  \item{methdesign}{
    Methylation study design matrix describing samples and groups. \strong{Use of edgeR::modelMatrixMeth() to make this matrix is highly recommended}, since it transforms a regular model.matrix (as one would construct for a microarray or RNA-Seq experiment) into a ``two-channel'' matrix representing methylated and unmethylated reads for each sample.
    Only applicable when \code{obj} is a BSseq object.
  }
  \item{all.cov}{
    If \code{TRUE}, only CpG sites where all samples have > 0 coverage will be retained. If \code{FALSE}, CpG sites for which some (not all) samples have coverage=0 will be retained. 
  }
  \item{contrasts}{
    Logical denoting whether a \code{limma}-style contrast matrix is specified.
    Only applicable when \code{obj} is a BSseq object.
  }
  \item{cont.matrix}{
    \code{Limma}-style contrast matrix for explicit contrasting. For each call to \code{sequencing.annotate}, only one contrast will be fit. 
    Only applicable when \code{obj} is a BSseq object.
  }
  \item{fdr}{
    FDR cutoff (Benjamini-Hochberg) for which CpG sites are individually called
    as significant. Used to index default thresholding in dmrcate(). \strong{Highly 
    recommended as the primary thresholding parameter for calling DMRs}. 
    Only applicable when \code{obj} is a BSseq object.
  }
  \item{coef}{
    The column index in \code{design} corresponding to the phenotype
    comparison. Corresponds to the comparison of interest in \code{design}
    when \code{contrasts=FALSE}, otherwise must be a column name in 
    \code{cont.matrix}. 
    Only applicable when \code{obj} is a BSseq object.
  }
  \item{\dots}{
    Extra arguments passed to the \code{limma} function lmFit().
    Only applicable when \code{obj} is a BSseq object.
  }
}
\value{
  A \code{\link{CpGannotated-class}}.
}
\references{
Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. \emph{Nucleic Acids Research}, \bold{43}(7), e47. 

Chen, Y., Pal, B., Visvader, J. E., & Smyth, G. K. (2017). Differential methylation analysis of reduced representation bisulfite sequencing experiments using edgeR. \emph{F1000Research}, \strong{6}, 2055. 

}

\author{
  Tim J. Peters <t.peters@garvan.org.au>
}
\examples{
library(ExperimentHub)
library(SummarizedExperiment)
eh = ExperimentHub()
bis_1072 <- eh[["EH1072"]]
pData(bis_1072) <- data.frame(replicate=gsub(".*-", "", colnames(bis_1072)),
                              tissue=substr(colnames(bis_1072), 1, nchar(colnames(bis_1072))-3), 
                              row.names=colnames(bis_1072))
colData(bis_1072)$tissue <- gsub("-", "_", colData(bis_1072)$tissue)
bis_1072 <- renameSeqlevels(bis_1072, mapSeqlevels(seqlevels(bis_1072), "UCSC"))
bis_1072 <- bis_1072[seqnames(bis_1072)=="chr19",]
bis_1072 <- bis_1072[240201:240300,]
tissue <- factor(pData(bis_1072)$tissue)
tissue <- relevel(tissue, "Liver_Treg")
design <- model.matrix(~tissue)
colnames(design) <- gsub("tissue", "", colnames(design))
colnames(design)[1] <- "Intercept"
rownames(design) <- colnames(bis_1072)
methdesign <- edgeR::modelMatrixMeth(design)
cont.mat <- limma::makeContrasts(treg_vs_tcon=Lymph_N_Treg-Lymph_N_Tcon,
                                 fat_vs_ln=Fat_Treg-Lymph_N_Treg,
                                 skin_vs_ln=Skin_Treg-Lymph_N_Treg,
                                 fat_vs_skin=Fat_Treg-Skin_Treg,
                                 levels=methdesign)
seq_annot <- sequencing.annotate(bis_1072, methdesign, all.cov = TRUE, 
                                   contrasts = TRUE, cont.matrix = cont.mat, 
                                   coef = "treg_vs_tcon", fdr=0.05)
}




