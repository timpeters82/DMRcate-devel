##
## get and set mixture-based ternary barcodes, usually at DMRs 
## see callDMRs for details and some cautions re: using a caller
## actually I should just set the TCGA AML DMR caller as the default...
## 
setGeneric('getBarcode', 
           function(object, ...) standardGeneric("getBarcode"))
setMethod('getBarcode', 'SummarizedExperiment', 
          function(object, ...) assays(object)$barcode)

setGeneric('setBarcode', 
           function(object, value, ...) standardGeneric("setBarcode"))
setMethod('setBarcode', 'SummarizedExperiment', 
          function(object, value, ...) assays(object)$barcode <- value)

