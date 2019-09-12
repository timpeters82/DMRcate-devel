\name{changeFDR}
\alias{changeFDR}
\title{
Change the individual CpG FDR thresholding for a CpGannotated object.
}
\description{
Takes a \code{\link{CpGannotated-class}} object and a specified FDR > 0 and < 1, and re-indexes the object in order to call DMRs at the specified rate.
}
\usage{
changeFDR(annot, FDR)
}
\arguments{
  \item{annot}{A \code{\link{CpGannotated-class}} object, created by \link{cpg.annotate} or \link{sequencing.annotate}.}
  \item{FDR}{The desired individual CpG FDR, which will index the rate at which DMRs are called.}
}
\details{
The number of CpG sites called as significant by this function will set the post-smoothing threshold for DMR constituents in \code{\link{dmrcate}}.
}
\value{
A re-indexed \code{\link{CpGannotated-class}} object.
}
\author{
Tim Peters <t.peters@garvan.org.au>
}
\examples{
data(dmrcatedata)
myMs <- logit2(myBetas)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
patient <- factor(sub("-.*", "", colnames(myMs)))
type <- factor(sub(".*-", "", colnames(myMs)))
design <- model.matrix(~patient + type) 
myannotation <- cpg.annotate("array", myMs.noSNPs, what="M", arraytype = "450K",
                             analysis.type="differential", design=design, coef=39)
myannotation_0.01 <- changeFDR(myannotation, 0.01)
}

