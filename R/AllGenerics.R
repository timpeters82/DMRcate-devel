setGeneric("changeFDR",  valueClass = "CpGannotated", function(annot, FDR) {
  standardGeneric("changeFDR")
})

setGeneric("cpg.annotate",  valueClass = "CpGannotated", function(datatype = c("array", "sequencing"), object, what = c("Beta", "M"), arraytype = c("EPIC", "450K"),
                                                                  analysis.type = c("differential", "variability", "ANOVA", "diffVar"), design, contrasts = FALSE, 
                                                                  cont.matrix = NULL, fdr = 0.05, coef, ...) {
  standardGeneric("cpg.annotate")
})

setGeneric("DMR.plot", function(ranges, dmr, CpGs, what = c("Beta", "M"), arraytype = c("EPIC", "450K"), phen.col, genome = c("hg19", "hg38", "mm10"), ...) {
  standardGeneric("DMR.plot")
})

setGeneric("dmrcate", valueClass="DMResults", function(object, lambda = 1000, C = NULL, pcutoff = "fdr", 
                               consec = FALSE, conseclambda = 10, min.cpgs = 2) {
  standardGeneric("dmrcate")
})

setGeneric("extractRanges", valueClass="GRanges", function(dmrcoutput, genome=c("hg19", "hg38", "mm10")) {
  standardGeneric("extractRanges")
})

setGeneric("rmSNPandCH", function(object, dist=2, mafcut=0.05, and=TRUE, rmcrosshyb=TRUE, rmXY=FALSE) {
  standardGeneric("rmSNPandCH")
})

setGeneric("sequencing.annotate", valueClass = "CpGannotated", function(obj, methdesign, all.cov=FALSE, contrasts = FALSE, 
                                           cont.matrix = NULL, fdr = 0.05, coef, ...) {
  standardGeneric("sequencing.annotate")
})
