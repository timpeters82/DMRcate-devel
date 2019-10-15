setMethod("show", "CpGannotated",
          function (object) cat(paste0("CpGannotated object describing ",
                                       length(object@ranges), " CpG sites, with independent\nCpG threshold indexed at fdr=",
                                       round(max(object@ranges$ind.fdr[object@ranges$is.sig]), 2), " and ", 
                                       sum(object@ranges$is.sig), " significant CpG sites.\n"))
)

setMethod("show", "DMResults",
          function (object) cat(paste0("DMResults object with ", length(object@coord), " DMRs.\nUse extractRanges() to produce a GRanges object of these.\n"))
)

