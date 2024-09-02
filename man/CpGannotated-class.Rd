\name{CpGannotated-class}
\alias{CpGannotated-class}
\docType{class}
\title{An object summarising individual CpG sites fitted to a given model}
\description{
  An S4 class that stores output from either \code{\link{cpg.annotate}} or \code{\link{sequencing.annotate}}.
}
\section{Slots}{
  \describe{
    \item{\code{ranges}:}{A GRanges object, containing CpG-level information to be passed to \code{\link{dmrcate}}. Mcols of this object include:
    
    - stat: Per-CpG test statistic; \emph{t} if from \code{limma} or Wald if from \code{DSS} if using differential mode. Variance if using variability mode, sqrt(\emph{F}) if using ANOVA mode, \emph{t} if using diffVar mode.
    
    - rawpval: Raw \emph{p}-value from DMP fitting. 
    
    - diff: Methylation difference/coefficient. In beta space for \code{\link{cpg.annotate}} output and output passed from \code{DSS::DMLtest()}. In logit space for when a \code{BSseq} object is passed from \code{\link{sequencing.annotate}}. Not available for output passed from \code{DSS::DMLtest.multiFactor()}. Not applicable in variability, ANOVA or diffVar modes.
    
    - ind.fdr: False discovery rate as calculated on individual CpG sites.
    
    - is.sig: Logical determining whether a CpG site is individually significant or not. Can be adjusted using \code{\link{changeFDR}}.}
    \item{\code{betas}:}{A matrix of per-CpG beta values matching the annotated loci.}
    }
    
}
\section{Methods}{
    \code{CpGannotated} objects have a show method that describes the data therein.
}
\author{ 
Tim Peters <t.peters@garvan.org.au>
}
\keyword{classes}

