\name{cpg.annotate}
\alias{cpg.annotate}
\title{
  Annotate Illumina CpGs with their chromosome position and test statistic
}
\description{
Annotate a matrix/GenomicRatioSet representing EPICv2, EPICv1 or 450K data with probe weights and chromosomal position. Provides replicate filtering and remapping functions for EPICv2 probes.

}
\usage{

cpg.annotate(datatype = c("array", "sequencing"), object,  
             what = c("Beta", "M"), arraytype = c("EPICv2", "EPICv1", "EPIC", 
             "450K"), epicv2Remap = TRUE, analysis.type = c("differential", 
             "variability", "ANOVA", "diffVar"), design, contrasts = FALSE, 
             cont.matrix = NULL, fdr = 0.05, coef, varFitcoef = NULL, 
             topVarcoef = NULL, ...) 
}
\arguments{
  \item{datatype}{
    Character string representing the type of data being analysed.
  }
  \item{object}{
    Either: 
    
    - A matrix of \emph{M}-values, with unique Illumina probe IDs as 
    rownames and unique sample IDs as column names or,
    
    - A GenomicRatioSet, appropriately annotated.
  }
  \item{what}{
    Does the data matrix contain Beta or M-values? Not needed
    if object is a GenomicRatioSet.
  }
  \item{arraytype}{
    Is the data matrix sourced from EPIC or 450K data? Not needed
    if object is a GenomicRatioSet.
  }
  \item{epicv2Remap}{
    Logical indicating whether to remap 11,878 cross-hybridising EPICv2 probes 
    to their more likely CpG target (see Peters \emph{et al.} 2024).
  }
  
  \item{analysis.type}{
    \code{"differential"} for \code{dmrcate()} to return DMRs;
    \code{"variability"} to return VMRs;
    \code{"ANOVA"} to return "whole experiment" DMRs, incorporating
    all possible contrasts from the design matrix using the moderated
    \emph{F}-statistics;
    \code{"diffVar"} to return differentially variable methylated regions, 
    using the \code{missMethyl} package to generate \emph{t}-statistics.
    
  }
  \item{design}{
    Study design matrix. Identical context to differential analysis
    pipeline in \code{limma}. Must have an intercept if \code{contrasts=FALSE}.
    Applies only when 
    \code{analysis.type \%in\% c("differential", "ANOVA", "diffVar")}.
  }
  \item{contrasts}{
    Logical denoting whether a \code{limma}-style contrast matrix is specified.
    Only applicable when \code{datatype="array"} and \code{analysis.type \%in\% 
    c("differential", "diffVar")}.
  }
  \item{cont.matrix}{
    \code{Limma}-style contrast matrix for explicit contrasting. For each call to \code{cpg.annotate}, only one contrast will be fit. 
    Only applicable when \code{datatype="array"} and \code{analysis.type \%in\% c("differential", "diffVar")}.
  }
  \item{fdr}{
    FDR cutoff (Benjamini-Hochberg) for which CpG sites are individually called
    as significant. Used to index default thresholding in dmrcate(). \strong{Highly 
    recommended as the primary thresholding parameter for calling DMRs}. 
    Not used when \code{analysis.type == "variability"}.
  }
  \item{coef}{
    The column index in \code{design} corresponding to the phenotype
    comparison. Corresponds to the comparison of interest in \code{design}
    when \code{contrasts=FALSE}, otherwise must be a column name in 
    \code{cont.matrix}. 
    Only applicable when \code{analysis.type == "differential"}.
  }
   \item{varFitcoef}{
    The columns of the design matrix containing the comparisons to test for
    differential variability. If left \code{NULL}, will test all columns.
    Identical context to \code{missMethyl::varFit()}. Only applicable when
    \code{analysis.type \%in\% "diffVar"}.
  }
  \item{topVarcoef}{
    Column number or column name specifying which coefficient of the linear
    model fit is of interest. It should be the same coefficient that 
    the differential variability testing was performed on. Default is last 
    column of fit object. Identical context to \code{missMethyl::topVar()}.
    Only applicable when \code{analysis.type \%in\% "diffVar"}.
  }
  \item{\dots}{
    Extra arguments passed to the \code{limma} function lmFit() (\code{analysis.type="differential"}).
  }
}
\value{
  A \code{\link{CpGannotated-class}}.
  }
\references{
Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. \emph{Nucleic Acids Research}, \bold{43}(7), e47. 

Phipson, B., & Oshlack, A. (2014). DiffVar: a new method for detecting differential variability with application to methylation in cancer and aging. \emph{Genome Biol}, \emph{15}(9), 465. 

Peters T.J., Buckley M.J., Statham, A., Pidsley R., Samaras K., Lord R.V., Clark S.J. and  Molloy P.L. \emph{De novo} identification of differentially methylated regions in the human genome. \emph{Epigenetics & Chromatin} 2015, \bold{8}:6,  doi:10.1186/1756-8935-8-6.

Peters, T.J., Meyer, B., Ryan, L., Achinger-Kawecka, J., Song, J., Campbell, E.M., Qu, W., Nair, S., Loi-Luu, P., Stricker, P., Lim, E., Stirzaker, C., Clark, S.J. and Pidsley, R. (2024). Characterisation and reproducibility of the HumanMethylationEPIC v2.0 BeadChip for DNA methylation profiling. \emph{BMC Genomics}, \bold{25}, 251. doi:10.1186/s12864-024-10027-5.
}

\author{
  Tim J. Peters <t.peters@garvan.org.au>
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
myannotation <- cpg.annotate("array", tcell, arraytype = "EPICv1",
                             analysis.type="differential", design=design, coef=2)
}




