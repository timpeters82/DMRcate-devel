\name{CpGannotated-class}
\alias{CpGannotated-class}
\docType{class}
\title{An object summarising individual CpG sites fitted to a given model}
\description{
  An S4 class that stores output from either \code{\link{cpg.annotate}} or \code{\link{sequencing.annotate}}.
}
\section{Slots}{
  This class has seven slots, containing CpG-level information to be passed to \code{\link{dmrcate}}:
  \describe{
    \item{\code{ID}:}{ Either an Illumina probe ID or a unique ID generated from \code{\link{sequencing.annotate}}.}
    \item{\code{stat}:}{ Per-CpG test statistic; \emph{t} if from \code{limma} or Wald if from \code{DSS}.}
    \item{\code{CHR}:}{ Chromosome of CpG site.}
    \item{\code{pos}:}{ Coordinate of CpG site.}
    \item{\code{diff}:}{ Methylation difference/coefficient. In beta space for \code{\link{cpg.annotate}} output and output passed from \code{DSS::DMLtest()}. In logit space for when a \code{BSseq} object is passed from \code{\link{sequencing.annotate}}. Not available for output passed from \code{DSS::DMLtest.multiFactor()}.}
    \item{\code{ind.fdr}:}{ False discovery rate as calculated on individual CpG sites.}
    \item{\code{is.sig}:}{ Logical determining whether a CpG site is individually significant or not. Can be adjusted using \code{\link{changeFDR}}.}
    }
}
\section{Methods}{
    \code{CpGannotate} objects have a \code{show} method that describes the data therein.
}
\author{ 
Tim Peters <t.peters@garvan.org.au>
}
\keyword{classes}

