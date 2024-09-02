\name{getCollapsedBetas}
\alias{getCollapsedBetas}
\title{
Extract a beta matrix from a CpGannotated object.
}
\description{
Retrieves the matrix of per-CpG beta values matching the annotated loci. For EPICv2 data, these betas have been (optionally) remapped to more likely target, and collapsed to one value per CpG. See \code{\link{cpg.annotate}} for details on EPICv2 replicate filtering.
}
\usage{
getCollapsedBetas(annot, ranges=NULL)
}
\arguments{
  \item{annot}{A \code{\link{CpGannotated-class}} object, created by \link{cpg.annotate}.}
  \item{ranges}{A GRanges object, over which the desired beta values are subsetted.}
}

\details{
Only applicable to \code{CpGannotated} objects that have been created by \code{cpg.annotate}; for those created by \code{sequencing.annotate} please use \code{bsseq::getMeth}.
}
\value{
A matrix of beta values, with individual CpG coordinates as rownames.
}
\author{
Tim Peters <t.peters@garvan.org.au>
}
\examples{
library(AnnotationHub)
ah <- AnnotationHub()
EPICv2manifest <- ah[["AH116484"]]
object <- minfi::logit2(matrix(rbeta(10000, 3, 1), 1000, 10))
rownames(object) <- sample(rownames(EPICv2manifest), 1000)
type <- rep(c("Ctrl", "Treat"), each=5)
design <- model.matrix(~type)
myannotation <- cpg.annotate("array", object, what = "M", arraytype = "EPICv2",
                             analysis.type="differential", design=design, coef=2)
getCollapsedBetas(myannotation)
}

