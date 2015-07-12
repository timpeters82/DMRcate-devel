\name{makeBedgraphs}
\alias{makeBedgraphs}
\title{
  Outputs bedGraphs
}
\description{
  Makes bedGraphs, 1 per sample, each containing all significant regions
  found via \code{dmrcate}. Bedgraphs are written are written to the
  working directory.
}
\usage{
makeBedgraphs(dmrcoutput, betas,
    annotation=c(array="IlluminaHumanMethylation450k",
        annotation="ilmn12.hg19"),
    samps=NULL)
}
\arguments{
  \item{dmrcoutput}{
     An object of class \code{dmrcate.output}.
  }
  \item{betas}{
    Matrix of beta values to be converted to bedGraph rows, with unique
    Illumina probe IDs as rownames.
  }
  \item{annotation}{
    A vector describing the type of annotation from which to derive
    bedgraph output. Identical context to \code{minfi},
    i.e. \code{annotation <- annotation(minfiobject)} where
    \code{minfiobject} is a \code{[Genomic](Methyl|Ratio)Set)}. Default
    (ilmn12.hg) is recommended.
  }
  \item{samps}{
    Vector of samples to be converted to bedGraph files, corresponding
    to columns of \code{betas}. Default is all samples plotted.
  }
}
\value{
  Writes zero or more bedGraph files to the working directory. 
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
    design=design, coef=39)
dmrcoutput <- dmrcate(myannotation, lambda=1000)
makeBedgraphs(dmrcoutput=dmrcoutput, betas=myBetas, samps=c(1,39))
}
}

