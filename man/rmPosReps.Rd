\name{rmPosReps}
\alias{rmPosReps}
\title{Filter out position replicates from an EPICv2 beta- or M-matrix}
\description{Given a beta-value or M-value matrix with EPICv2 probe IDs as rownames, returns a truncated matrix with a 1-to-1 mapping of probe ID to CpG locus. Values returned depend on the filtering strategy selected. Replicate probes are averaged by default, but the user may optionally select individual probes per replicate group based on maximum sensitivity to methylation change or maximum precision, as per empirical cross-platform consensus testing against EPICv1 and WGBS data (Peters \emph{et al.} 2024).}
\usage{
rmPosReps(object, filter.strategy= c("mean", "sensitivity",
                                     "precision","random"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
   A matrix of beta- or \emph{M}-values, with unique EPICv2 Illumina probe IDs as 
    rownames.
  }

  \item{filter.strategy}{
    Strategy for filtering probe replicates that map to the same CpG site.
    \code{"mean"} takes the mean of the available probes;
    \code{"sensitivity"} takes the available probe most sensitive to methylation
    change;
    \code{"precision"} either selects the available probe with the lowest 
    variation from the consensus value (most precise), or takes the mean if that
    confers the lowest variation instead,
    \code{"random"} takes a single probe at random from each replicate group.
  }
}
\value{A truncated matrix with a 1-to-1 mapping of probe to CpG site. If the group
mean is taken, the first probe in the group by alphabetical sorting is returned as
the rowname.
}
\references{
Peters, T.J., Meyer, B., Ryan, L., Achinger-Kawecka, J., Song, J., 
  Campbell, E.M., Qu, W., Nair, S., Loi-Luu, P., Stricker, P., Lim, E., 
  Stirzaker, C., Clark, S.J. and Pidsley, R. (2024). Characterisation and 
  reproducibility of the HumanMethylationEPIC v2.0 BeadChip for DNA methylation 
  profiling. \emph{BMC Genomics} \bold{25}, 251. doi:10.1186/s12864-024-10027-5.
}
\author{
Tim J. Peters <t.peters@garvan.org.au>
}

\examples{
library(ExperimentHub)
eh <- ExperimentHub()
ALLbetas <- eh[["EH9451"]]
ALLbetas <- ALLbetas[1:1000,]
ALLMs <- minfi::logit2(ALLbetas)
ALLMs.repmean <- rmPosReps(ALLMs, filter.strategy="mean")
}

