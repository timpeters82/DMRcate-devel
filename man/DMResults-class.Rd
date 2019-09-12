\name{DMResults-class}
\alias{DMResults-class}
\docType{class}
\title{Initial storage object for called DMRs - class}
\description{
An S4 class that stores DMR information as output from \code{\link{dmrcate}}.
}
\section{Slots}{
This class has eight slots, summarising DMR information to be passed to \code{\link{extractRanges}}:
  \describe{
    \item{\code{coord}:}{ DMR coordinates in UCSC style.}
    \item{\code{no.cpgs}:}{ Number of constituent CpG sites of DMR.}
    \item{\code{min_smoothed_fdr}:}{ Minimum FDR of the \strong{smoothed} estimate.}
    \item{\code{Stouffer}:}{ Stouffer summary transform of the \strong{individual} CpG FDRs.}
    \item{\code{HMFDR}:}{ Harmonic mean of the \strong{individual} CpG FDRs.}
    \item{\code{Fisher}:}{ Fisher combined probability transform of the \strong{individual} CpG FDRs.}
    \item{\code{maxdiff}:}{ Maximum differential/coefficient within the DMR.}
    \item{\code{meandiff}:}{ Mean differential/coefficient across the DMR.}
    }
}
\section{Methods}{
\code{DMResults} objects have a \code{show} method describing the number of DMRs called.
}
\author{Tim Peters <t.peters@garvan.org.au>}
\keyword{classes}

