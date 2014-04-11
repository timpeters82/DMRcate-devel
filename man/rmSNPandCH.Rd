\name{rmSNPandCH}
\alias{rmSNPandCH}
\title{
  Filter probes
}
\description{
  Filters a matrix of M-values (or beta values) by distance.
}
\usage{
rmSNPandCH(object, dist = 2, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE)
}
\arguments{
  \item{object}{
    A matrix of M-values or beta values, with unique Illumina probe IDs
    as rownames.
  }
  \item{dist}{
    Maximum distance (from CpG to SNP) of probes to be filtered out. See
    details for when Illumina occasionally lists a CpG-to-SNP distance
    as being < 0.
  }
  \item{mafcut}{
    Minimum minor allele frequency of probes to be filtered out.
  }
  \item{and}{
    If \code{TRUE}, the probe must have at least 1 SNP binding to it
    that satisfies \bold{both} requirements in \code{dist} and
    \code{mafcut} for it to be filtered out. If \code{FALSE}, it will be
    filtered out if either requirement is satisfied. Default is
    \code{TRUE}.
  }
  \item{rmcrosshyb}{
    If \code{TRUE}, filters out probes found by Chen et al. (2013) to be
    cross-reactive with areas of the genome not at the site of
    interest. Many of these sites are on the X-chromosome, leading to
    potential confounding if the sample group is a mix of males and
    females. There are 30,969 probes in total in this list. Default is
    \code{TRUE}.
  }
}
\details{
  Probes in \code{-1:dist} will be filtered out for any integer
  specification of \code{dist}. When a probe is listed as being ``-1''
  nucleotides from a SNP (7 in total of the 153,113), that SNP is
  immediately adjacent to the end of the probe, and is likely to
  confound the measurement, in addition to those listed as 0, 1 or 2
  nucleotides away. See vignette for further details.
}
\value{
  A matrix, attenuated from \code{object}, with rows corresponding to
  probes matching user input filtered out.
}
\references{
  Chen YA, Lemire M, Choufani S, Butcher DT, Grafodatskaya D, Zanke BW,
  Gallinger S, Hudson TJ, Weksberg R. Discovery of cross-reactive probes
  and polymorphic CpGs in the Illumina Infinium HumanMethylation450
  microarray. \emph{Epigenetics}. 2013 Jan 11;8(2).

  \url{http://supportres.illumina.com/documents/myillumina/88bab663-307c-444a-848e-0ed6c338ee4d/humanmethylation450_15017482_v.1.2.snpupdate.table.v3.txt}

}
\author{
  Tim J. Peters <Tim.Peters@csiro.au>
}
\examples{
\dontrun{
data(dmrcatedata)
myMs <- logit2(myBetas)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
}
}
