\name{extractRanges}
\alias{extractRanges}
\title{
Create GRanges object from \code{\link{dmrcate}} output.
}
\description{
Takes a \code{dmrcate.output} object and produces the corresponding GRanges object. 
}
\usage{
extractRanges(dmrcoutput)
}

\arguments{
  \item{dmrcoutput}{
    An object of class \code{dmrcate.output}.
}
}
\value{
A GRanges object.
}
\author{
Tim Triche Jr. <tim.triche@usc.edu>, Tim Peters <Tim.Peters@csiro.au>
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
myRanges <- extractRanges(dmrcoutput)
}
}
