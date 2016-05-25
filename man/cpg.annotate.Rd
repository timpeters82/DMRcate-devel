\name{cpg.annotate}
\alias{cpg.annotate}
\title{
  Annotate CpGs with their chromosome position and test statistic
}
\description{
  Either:
  - Annotate a matrix of \emph{M}-values (logit transform of beta)
  representing 450K or EPIC data with probe weights (depending on \code{analysis.type}) 
  and chromosomal position, or
  - Standardise this information from \code{DSS:::DMLtest()} to the
    same data format.
  
}
\usage{

cpg.annotate(datatype = c("array", "sequencing"), 
             object, 
             annotation=c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"), 
             analysis.type = c("differential", "variability"), 
             design, 
             contrasts = FALSE, 
             cont.matrix = NULL, 
             fdr = 0.05, 
             coef, 
             ...) 
}
\arguments{
  \item{datatype}{
    Character string representing the type of data being analysed.
  }
  \item{object}{
    Either: 
    
    - A matrix of \emph{M}-values, with unique Illumina probe IDs as 
    rownames and unique sample IDs as column names or, 
    
    - Output from \code{DSS:::DMLtest()}.
  }
  \item{annotation}{
    A vector describing the type of annotation to affix to
    \code{object}. Identical context to \code{minfi},
    i.e. \code{annotation <- annotation(minfiobject)} where
    \code{minfiobject} is a \code{[Genomic](Methyl|Ratio)Set)}. 
    
    Argument for 450K arrays: 
    
    \code{c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19")}.
    
    Argument for EPIC arrays:
    
    \code{c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")}.
    
    An error will be thrown if you attempt one on an \code{object}
    with rownames on the other, due to non-overlapping probes 
    on both platforms. Only applicable when \code{datatype="array"}.
  }
  \item{analysis.type}{
    \code{"differential"} for \code{dmrcate()} to return DMRs and
    \code{"variability"} to return VMRs. Only applicable when \code{datatype="array"}.
  }
  \item{design}{
    Study design matrix. Identical context to differential analysis
    pipeline in \code{limma}. Must have an intercept if \code{contrasts=FALSE}.
    Applies only when \code{analysis.type="differential"}. 
    Only applicable when \code{datatype="array"}.
  }
  \item{contrasts}{
    Logical denoting whether a \code{limma}-style contrast matrix is specified.
    Only applicable when \code{datatype="array"}.
  }
  \item{cont.matrix}{
    \code{Limma}-style contrast matrix for explicit contrasting. For each call to \code{cpg.annotate}, only one contrast will be fit. 
    Only applicable when \code{datatype="array"}.
  }
  \item{fdr}{
    FDR cutoff (Benjamini-Hochberg) for which CpG sites are individually called
    as significant. Used only to determine effect size, and not for
    downstream thresholding.
  }
  \item{coef}{
    The column index in \code{design} corresponding to the phenotype
    comparison. Corresponds to the comparison of interest in \code{design}
    when \code{contrasts=FALSE}, otherwise must be a column name in 
    \code{cont.matrix}. Applies only when \code{analysis.type="differential"}
    and when \code{datatype="array"}.
  }
  \item{\dots}{
    Extra arguments passed to the \code{limma} function lmFit(). \cr
    Applies only when \code{analysis.type="differential"} and when \code{datatype="array"}.
  }
}
\value{
  An object of class "annot", for passing to \code{dmrcate}, containing
  the vectors:
  \itemize{
    \item \code{ID}: Illumina probe ID or row number
    \item \code{stat}: \emph{t}-statistic or Wald statistics between phenotypes for each CpG
    \item \code{CHR}: Chromosome which the CpG maps to
    \item \code{pos}: Genomic coordinate (on \code{CHR}) that the CpG maps to
    \item \code{betafc}: The beta fold change according to the given design
    \item \code{indfdr}: Individually-derived FDRs for each CpG 
    \item \code{is.sig}: Logical denoting either significance from \code{fdr} (\code{analysis.type="differential"}) or top ventile of variable probes (\code{analysis.type="variability"})
  }
}
\references{
Smyth, G. K. (2005). Limma: linear models for microarray data. In: \emph{Bioinformatics and Computational Biology Solutions using R and Bioconductor}, R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds.), Springer, New York, pages 397-420.

Feng, H., Conneely, K. N., & Wu, H. (2014). A Bayesian hierarchical model to detect differentially methylated loci from single nucleotide resolution sequencing data. \emph{Nucleic Acids Research}, \bold{42}(8), e69. 

Peters T.J., Buckley M.J., Statham, A., Pidsley R., Samaras K., Lord R.V., Clark S.J. and  Molloy P.L. \emph{De novo} identification of differentially methylated regions in the human genome. \emph{Epigenetics & Chromatin} 2015, \bold{8}:6,  doi:10.1186/1756-8935-8-6.
}

\author{
  Tim J. Peters <t.peters@garvan.org.au>
}
\examples{
\dontrun{
data(dmrcatedata)
myMs <- logit2(myBetas)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
patient <- factor(sub("-.*", "", colnames(myMs)))
type <- factor(sub(".*-", "", colnames(myMs)))
design <- model.matrix(~patient + type) 
myannotation <- cpg.annotate("array", myMs.noSNPs, analysis.type="differential",
    design=design, coef=39)
}
}




