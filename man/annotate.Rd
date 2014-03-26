\name{annotate}
\alias{annotate}
\title{
  450k probe annotation
}
\description{
  Annotates a matrix of \emph{M}-values (logit transform of beta) with weights
  (depending on \code{analysis.type}) and other relevant information
  including gene association.
}
\usage{

annotate(object,
    annotation=c(array="IlluminaHumanMethylation450k",
        annotation="ilmn12.hg19"),
    analysis.type=c("differential", "variability", "hypermethylation"),
    design, coef, diff.metric=c("FC", "t"), paired=FALSE, pcutoff=NULL,
    betacutoff=NULL, quantcut=0.9, cut=0.5, ...)
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
    \code{"differential"} for \code{dmrcate()} to return DMRs,
    \code{"variability"} to return VMRs and \code{"hypermethylation"} to
    return HMRs.
  }
  \item{design}{
    Study design matrix. Identical context to differential analysis
    pipeline in \code{limma}. Applies only when
    \code{analysis.type="differential"}.
  }
  \item{coef}{
    The column index in \code{design} corresponding to the phenotype
    comparison. Identical context to the \code{coef} argument in
    \code{limma}'s \code{topTable}. Applies only when
    \code{analysis.type="differential"}.
  }
  \item{diff.metric}{
    The value with which the density estimate is weighted with in
    \code{dmrcate}: either fold-change (\code{"FC"}) or
    \emph{t}-statistic (\code{"t"}).  Applies only when
    \code{analysis.type="differential"}.
  }

  \item{paired}{
    Logical; denotes whether study design has paired samples
    (i.e. biological replicates) or not. This should be reflected in
    \code{design}. \cr Applies only when
    \code{analysis.type="differential"}.
  }
  \item{pcutoff}{
    If specified, the BH-adjusted \emph{p}-value threshold at which probes are
    filtered for passing to \code{dmrcate}. Must be a numeric argument 0
    <= \code{pcutoff} <= 1. Default is \code{NULL} i.e. no
    filtering. User should be confident that enough probes will survive
    the significance test so that \code{dmrcate} can fit a meaningful
    density estimate. If in doubt, the default should be used. Applies
    only when \code{analysis.type="differential"}.
  }
  \item{betacutoff}{
    If specified, the minimum beta fold change at which probes are
    filtered for passing to \code{dmrcate}. Cannot be used if
    \code{pcutoff} is specified. The same caveats for specifying
    \code{pcutoff} also apply here. \cr Applies only when
    \code{analysis.type="differential"}.
  }
  \item{quantcut}{
    Proportion of probes to be filtered out by variance. \cr Applies
    only when \code{analysis.type="variability"}.
  }
  \item{cut}{
    Cutoff for retaining only hypermethylated probes. Recommended
    default is 0.5. Applies only when
    \code{analysis.type="hypermethylation"}.
  }
  \item{\dots}{
    Extra arguments passed to the \code{limma} function lmFit(). \cr
    Applies only when \code{analysis.type="differential"}.
  }
}
\details{
  Specification of a \code{pcutoff} value often results in very
  conservative estimates for significance of DMRs when the resulting
  \code{DMannot} object is passed to \code{dmrcate}. When the two
  phenotypes being compared are obviously very different (for example
  cancer vs. matched normal) then it is probably safe to make
  \code{pcutoff=0.05}. However, for subtle phenotype comparisons,
  especially those that are non-matched, we recommend using the default.
}
\value{
  An object of class "annot", for passing to \code{dmrcate}, containing
  the vectors:
  \itemize{
    \item \code{ID}: Illumina probe ID
    \item \code{weights}: Fold change between phenotypes for each probe
    \item \code{CHR}: Chromosome which the probe maps to
    \item \code{pos}: hg19 position (on \code{CHR}) that the probe maps to
    \item \code{gene}: Matching UCSC_RefGene_Name
    \item \code{group}: Matching UCSC_RefGene_Group
    \item \code{betafc}: The beta fold change according to the given design
  }
}
\references{
  Smyth, G. K. (2004). Linear models and empirical Bayes methods for
  assessing differential expression in microarray
  experiments. \emph{Statistical Applications in Genetics and Molecular
  Biology} \bold{3}, No. 1, Article 3.
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
    design=design, coef=39, diff.metric="FC", paired=TRUE, pcutoff=0.01)}
}




