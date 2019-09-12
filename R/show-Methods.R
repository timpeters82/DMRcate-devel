setMethod("show", "CpGannotated",
          function (object) cat(paste0("CpGannotated object describing ",
                                       length(object@ID), " CpG sites, with independent\nCpG threshold indexed at fdr=",
                                       round(max(object@ind.fdr[object@is.sig]), 2), " and ", 
                                       sum(object@is.sig), " significant CpG sites.\n"))
)

setMethod("show", "DMResults",
          function (object) cat(paste0("DMResults object with ", length(object@coord), " DMRs.\nUse extractRanges() to produce a GRanges object of these.\n"))
)

