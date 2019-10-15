setClass("CpGannotated",
         slots = c(ranges="GRanges"),
         prototype = list(ranges=GRanges()))
setClass("DMResults",
         slots = c(coord="character", no.cpgs="integer", min_smoothed_fdr="numeric", Stouffer="numeric", 
                   HMFDR="numeric", Fisher="numeric", maxdiff="numeric", meandiff="numeric"),
         prototype = list(coord=character(), no.cpgs=integer(), min_smoothed_fdr=numeric(), Stouffer=numeric(),
                   HMFDR=numeric(), Fisher=numeric(), maxdiff=numeric(), meandiff=numeric()))

