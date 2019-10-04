\name{extractRanges}
\alias{extractRanges}
\title{
Create a GRanges object from \code{\link{dmrcate}} output.
}
\description{
Takes a DMResults object and produces the corresponding GRanges object. 
}
\usage{
extractRanges(dmrcoutput, genome = c("hg19", "hg38", "mm10")) 
}

\arguments{
  \item{dmrcoutput}{
    A DMResults object.
}
  \item{genome}{
    Reference genome for annotating DMRs with promoter overlaps. 
    Can be one of \code{"hg19"}, \code{"hg38"} or \code{"mm10"}.
    Ranges are assumed to map to the reference stated; there is no liftover.
  }
}
\value{
A GRanges object.
}
\author{
Tim Triche Jr. <tim.triche@usc.edu>, Tim Peters <t.peters@garvan.org.au>
}
\examples{
library(ExperimentHub)
eh <- ExperimentHub()
FlowSorted.Blood.EPIC <- eh[["EH1136"]]  
tcell <- FlowSorted.Blood.EPIC[,colData(FlowSorted.Blood.EPIC)$CD4T==100 | 
                                 colData(FlowSorted.Blood.EPIC)$CD8T==100]
locs <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
locs <- locs[order(locs$chr, locs$pos),]
chr2probes <- rownames(locs)[457851:457950]
tcell <- subsetByLoci(tcell, includeLoci = chr2probes)
detP <- detectionP(tcell)
remove <- apply(detP, 1, function (x) any(x > 0.01))
tcell <- tcell[!remove,]
tcell <- preprocessNoob(tcell)
tcellms <- getM(tcell)
tcellms.noSNPs <- rmSNPandCH(tcellms, dist=2, mafcut=0.05)
tcell$Replicate[tcell$Replicate==""] <- tcell$Sample_Name[tcell$Replicate==""]
tcellms.noSNPs <- limma::avearrays(tcellms.noSNPs, tcell$Replicate)
tcell <- tcell[,!duplicated(tcell$Replicate)]
type <- factor(tcell$CellType)
design <- model.matrix(~type) 
myannotation <- cpg.annotate("array", tcellms.noSNPs, what="M", arraytype = "EPIC",
                             analysis.type="differential", design=design, coef=2)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
}
