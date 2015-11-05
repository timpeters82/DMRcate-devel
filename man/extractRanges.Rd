\name{extractRanges}
\alias{extractRanges}
\title{
Create GRanges object from \code{\link{dmrcate}} output.
}
\description{
Takes a \code{dmrcate.output} object and produces the corresponding GRanges object. 
}
\usage{
extractRanges(dmrcoutput, genome = c("hg19", "hg38", "mm10")) 
}

\arguments{
  \item{dmrcoutput}{
    An object of class \code{dmrcate.output}.
}
  \item{genome}{
    Reference genome for annotating DMRs with promoter overlaps. 
    Can be one of \code{"hg19"}, \code{"hg38"} or \code{"mm10"}
  }
}
\value{
A GRanges object.
}
\author{
Tim Triche Jr. <tim.triche@usc.edu>, Tim Peters <t.peters@garvan.org.au>
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
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
}
}
