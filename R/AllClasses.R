setClass("CpGannotated",
         slots = c(ID="character", stat="numeric", CHR="character", pos="integer", 
                   diff="numeric", ind.fdr="numeric", is.sig="logical"),
         prototype = list(ID=character(), stat=numeric(), CHR=character(), pos=integer(),
                   diff=numeric(), ind.fdr=numeric(), is.sig=logical()))
setClass("DMResults",
         slots = c(coord="character", no.cpgs="integer", min_smoothed_fdr="numeric", Stouffer="numeric", 
                   HMFDR="numeric", Fisher="numeric", maxdiff="numeric", meandiff="numeric"),
         prototype = list(coord=character(), no.cpgs=integer(), min_smoothed_fdr=numeric(), Stouffer=numeric(),
                   HMFDR=numeric(), Fisher=numeric(), maxdiff=numeric(), meandiff=numeric()))

