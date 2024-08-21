rmSNPandCH <- function (object, dist = 2, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY = FALSE) 
{
  stopifnot(is.matrix(object))
  dist <- as.integer(dist)
  stopifnot(0 <= mafcut & mafcut <= 1)
  
  if(nchar(rownames(object))[1] < 10){
    stop("Error: rownames of object do not seem to be Illumina probe IDs. Please annotate your object or assign probe IDs to rownames.")
  }
  
  if(nchar(rownames(object))[1] == 10){
    
    message("Probe IDs from EPICv1 or earlier detected. Proceeding...")
    eh <- ExperimentHub()
    snpsall <- eh[["EH3130"]]
    dist0 <- as.integer(snpsall$distances)
    distrange <- range(dist0)
    stopifnot(dist >= min(distrange) && dist <= max(distrange))
    test0 <- (dist0 >= -1) & (dist0 <= dist)
    test1 <- snpsall$mafs > mafcut
    test <- if (and) 
      (test0 & test1)
    else (test0 | test1)
    badprobes <- snpsall$probe[test]
    if (rmcrosshyb) {
      crosshyb <- eh[["EH3129"]]
      badprobes <- union(badprobes, as.character(crosshyb))
    }
    if (rmXY) {
      XY.probes <- eh[["EH3131"]]
      badprobes <- union(badprobes, XY.probes)
    }
  } else {
    message("Assuming EPICv2 data. Proceeding...")
    eh <- ExperimentHub()
    epicv2snps <- eh[["EH8568"]]
    dist0 <- as.integer(epicv2snps$distances)
    distrange <- range(dist0)
    stopifnot(dist >= min(distrange) && dist <= max(distrange))
    test0 <- (dist0 >= 0) & (dist0 <= dist)
    test1 <- epicv2snps$mafs > mafcut
    test <- if (and) { 
      (test0 & test1)
    }  else (test0 | test1) 
    badprobes <- epicv2snps$probe[test]
    if (rmcrosshyb) {
      ah <- AnnotationHub()
      EPICv2manifest <- ah[["AH116484"]]
      crosshyb <- rownames(EPICv2manifest)[EPICv2manifest$CH_BLAT=="Y"]
      badprobes <- union(badprobes, as.character(crosshyb))
    }
    if (rmXY) {
      if(!exists("EPICv2manifest")){
        ah <- AnnotationHub()
        EPICv2manifest <- ah[["AH116484"]]
      }
      XY.probes <- rownames(EPICv2manifest)[EPICv2manifest$CHR %in% c("chrX", "chrY")]
      badprobes <- union(badprobes, XY.probes)
    }
  }
  #Only retain cytosines
  keep <- grep("^cg|^ch", rownames(object))
  object <- object[keep,]
  
  object[!(rownames(object) %in% badprobes), ]
}
