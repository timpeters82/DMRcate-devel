\name{rmSNPandCH}
\alias{rmSNPandCH}
\title{
  Filter probes
}
\description{
  Filters a matrix of M-values (or beta values) by distance to SNP/variant. Also (optionally) removes cross-hybridising probes and sex-chromosome probes.
}
\usage{
rmSNPandCH(object, dist = 2, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY=FALSE)
}
\arguments{
  \item{object}{
    A matrix of M-values or beta values, with unique Illumina probe IDs
    as rownames.
  }
  \item{dist}{
    Maximum distance (from CpG to SNP/variant) of probes to be filtered out. See
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
    If \code{TRUE}, filters out probes found by Pidsley and Zotenko et al. (2016)
    for EPIC or Chen et al. (2013) for 450K to be cross-reactive with areas of the genome not
    at the site of interest. Many of these sites are on the X-chromosome, leading to
    potential confounding if the sample group is a mix of males and
    females. There are 63,707 probes in total in this list. Default is
    \code{TRUE}.
  }
  \item{rmXY}{
    If \code{TRUE}, filters out probe hybridising to sex chromosomes. Or-operator applies when combined with other 2 filters. 
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
  Pidsley R, Zotenko E, Peters TJ, Lawrence MG, Risbridger GP, Molloy P, 
  Van Dijk S, Muhlhausler B, Stirzaker C, Clark SJ. Critical evaluation of 
  the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA 
  methylation profiling. \emph{Genome Biology}. 2016 17(1), 208.
  
  Chen YA, Lemire M, Choufani S, Butcher DT, Grafodatskaya D, Zanke BW,
  Gallinger S, Hudson TJ, Weksberg R. Discovery of cross-reactive probes
  and polymorphic CpGs in the Illumina Infinium HumanMethylation450
  microarray. \emph{Epigenetics}. 2013 Jan 11;8(2).
  \url{http://supportres.illumina.com/documents/myillumina/88bab663-307c-444a-848e-0ed6c338ee4d/humanmethylation450_15017482_v.1.2.snpupdate.table.v3.txt}

}
\author{
  Tim J. Peters <t.peters@garvan.org.au>
}
\examples{
library(ExperimentHub)
eh <- ExperimentHub()
FlowSorted.Blood.EPIC <- eh[["EH1136"]]  
tcell <- FlowSorted.Blood.EPIC[,colData(FlowSorted.Blood.EPIC)$CD4T==100 | 
                                 colData(FlowSorted.Blood.EPIC)$CD8T==100]
                                 detP <- detectionP(tcell)
detP <- detectionP(tcell)
remove <- apply(detP, 1, function (x) any(x > 0.01))
tcell <- tcell[!remove,]
tcell <- preprocessFunnorm(tcell)
tcellms <- getM(tcell)
tcellms.noSNPs <- rmSNPandCH(tcellms, dist=2, mafcut=0.05)

}
