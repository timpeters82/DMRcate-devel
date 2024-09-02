setClass("CpGannotated",
         slots = c(ranges="GRanges", betas="matrix"),
         prototype = list(ranges=GRanges(), betas=matrix()))
setClass("DMResults",
         slots = c(coord="character", no.cpgs="integer", min_smoothed_fdr="numeric", Stouffer="numeric", 
                   HMFDR="numeric", Fisher="numeric", maxdiff="numeric", meandiff="numeric"),
         prototype = list(coord=character(), no.cpgs=integer(), min_smoothed_fdr=numeric(), Stouffer=numeric(),
                   HMFDR=numeric(), Fisher=numeric(), maxdiff=numeric(), meandiff=numeric()))

