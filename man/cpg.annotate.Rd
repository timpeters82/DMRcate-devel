\name{cpg.annotate}
\alias{cpg.annotate}
\title{
  450k probe annotation
}
\description{
  Annotates a matrix of \emph{M}-values (logit transform of beta) with weights
  (depending on \code{analysis.type}) and other relevant information
  including gene association.
}
\usage{

cpg.annotate(object, annotation=c(array="IlluminaHumanMethylation450k",
                                      annotation="ilmn12.hg19"),
         analysis.type=c("differential", "variability"),
         design, contrasts=FALSE, cont.matrix=NULL, coef, ...)
}
\arguments{
  \item{object}{
    A matrix of M-values, with unique Illumina probe IDs as rownames and
    unique sample IDs as column names.
  }
  \item{annotation}{
    A vector describing the type of annotation to affix to
    \code{object}. Identical context to \code{minfi},
    i.e. \code{annotation <- annotation(minfiobject)} where
    \code{minfiobject} is a \code{[Genomic](Methyl|Ratio)Set)}. Default
    (ilmn12.hg) is recommended.
  }
  \item{analysis.type}{
    \code{"differential"} for \code{dmrcate()} to return DMRs and
    \code{"variability"} to return VMRs.
  }
  \item{design}{
    Study design matrix. Identical context to differential analysis
    pipeline in \code{limma}. Must have an intercept if \code{contrasts=FALSE}.
    Applies only when \code{analysis.type="differential"}.
  }
  \item{contrasts}{
    Logical denoting whether a \code{limma}-style contrast matrix is specified.
  }
  \item{cont.matrix}{
    \code{Limma}-style contrast matrix for explicit contrasting. For each call to \code{cpg.annotate}, only one contrast will be fit. 
  }
  \item{coef}{
    The column index in \code{design} corresponding to the phenotype
    comparison. Corresponds to the comparison of interest in \code{design}
    when \code{contrasts=FALSE}, otherwise must be a column name in 
    \code{cont.matrix}. Applies only when \code{analysis.type="differential"}.
  }
  \item{\dots}{
    Extra arguments passed to the \code{limma} function lmFit(). \cr
    Applies only when \code{analysis.type="differential"}.
  }
}
\value{
  An object of class "annot", for passing to \code{dmrcate}, containing
  the vectors:
  \itemize{
    \item \code{ID}: Illumina probe ID
    \item \code{weights}: \emph{t}-statistic between phenotypes for each probe
    \item \code{CHR}: Chromosome which the probe maps to
    \item \code{pos}: hg19 position (on \code{CHR}) that the probe maps to
    \item \code{gene}: Matching UCSC_RefGene_Name
    \item \code{group}: Matching UCSC_RefGene_Group
    \item \code{betafc}: The beta fold change according to the given design
    \item \code{indfdr}: The post-kernel fitting \code{limma} fdr value 
  }
}
\references{
Smyth, G. K. (2005). Limma: linear models for microarray data. In: \emph{Bioinformatics and Computational Biology Solutions using R and Bioconductor}, R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds.), Springer, New York, pages 397-420.

Peters T.J., Buckley M.J., Statham, A., Pidsley R., Samaras K., Lord R.V., Clark S.J. and  Molloy P.L. \emph{De novo} identification of differentially methylated regions in the human genome. \emph{Epigenetics & Chromatin} 2015, \bold{8}:6,  doi:10.1186/1756-8935-8-6.
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
myannotation <- cpg.annotate(myMs.noSNPs, analysis.type="differential",
    design=design, coef=39)}
}




