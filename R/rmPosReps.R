rmPosReps <- function(object, filter.strategy= c("mean", "sensitivity","precision","random")){ 
  if(any(nchar(rownames(object)) < 13)){
    stop("Error: rownames do not look like EPICv2 probes. This function will only work for EPICv2 data.")
  }
  if(any(!grepl("^cg|^ch", rownames(object)))){
    stop("Error: This function will only accept a matrix with rownames beginning with cg or ch. Please run your matrix through rmSNPandCH() first.")
  }
  message("Loading EPICv2 manifest...")
  ah <- AnnotationHub()
  EPICv2manifest <- ah[["AH116484"]]
  EPICv2manifest <- EPICv2manifest[rownames(object),]
  EPICv2manifest <- EPICv2manifest[EPICv2manifest$CHR!="chr0",]
  object <- object[rownames(EPICv2manifest),]
  coords <- paste(EPICv2manifest$CHR, EPICv2manifest$MAPINFO, sep=":")
  posreps <- table(coords)
  posreps <- names(posreps)[posreps > 1]
  if(length(posreps)==0){
    stop("No EPICv2 replicates detected. Please proceed to cpg.annotate().")
  }
  
  switch(filter.strategy, mean={
    message("Averaging probes that map to the same CpG site...")
    outs <- lapply(posreps, function (x){
      ids <- coords==x
      means <- colMeans(object[ids,])
      retain <- rownames(EPICv2manifest)[ids][1]
      dups <- rownames(EPICv2manifest)[ids][-1]
      list(means, retain, dups)
    })
    means <- do.call("rbind", lapply(outs, function (x) x[[1]]))
    rownames(means) <- unlist(lapply(outs, function (x) x[[2]]))
    object[rownames(means),] <- means
    dups <- unlist(lapply(outs, function (x) x[[3]]))
    EPICv2manifest <- EPICv2manifest[!rownames(EPICv2manifest) %in% dups,]
    object <- object[!rownames(object) %in% dups,]
  }, sensitivity={
    message("Selecting probes that map to the same CpG site by sensitivity to methylation change...")
    senschoice <- lapply(posreps, function (x) {
      probes <- rownames(EPICv2manifest)[coords==x]
      classes <- EPICv2manifest[probes, "Rep_results_by_LOCATION"]
      if(any(grepl("Superior|sensitivity", classes))){
        if(length(grep("Superior|sensitivity", classes)) > 1){
          choice <- sample(probes[grepl("Superior|sensitivity", classes)], 1)
        } else {
          choice <- probes[grepl("Superior|sensitivity", classes)]
        }
      } 
      else if (any(grepl("Inferior", classes))){
        if(all(grepl("Inferior", classes))){
          choice <- sample(probes, 1)
        } else {
          choice <- sample(probes[!grepl("Inferior", classes)], 1)
        }
      } else {
        choice <- sample(probes, 1)
      }
      dups <- probes[!probes==choice]
      list(choice, dups)
    })
    remove <- unlist(lapply(senschoice, function (x) x[[2]]))
    EPICv2manifest <- EPICv2manifest[!rownames(EPICv2manifest) %in% remove,]
    object <- object[!rownames(object) %in% remove,]
  }, precision={
    message("Processing probes that map to the same CpG site for best precision...")
    precchoice <- lapply(posreps, function (x) {
      probes <- rownames(EPICv2manifest)[coords==x]
      classes <- EPICv2manifest[probes, "Rep_results_by_LOCATION"]
      if(any(grepl("mean", classes))){
        means <- colMeans(object[probes,])
        retain <- probes[1]
        dups <- probes[-1]
        return(list(means, retain, dups))
      } 
      else if (any(classes=="Best precision")){
        if(sum(classes=="Best precision") > 1){
          choice <- sample(probes[classes=="Best precision"], 1)
        } else {
          choice <- probes[grepl("Best precision", classes)]
        }
      } 
      else if (any(grepl("Superior", classes))){
        if(length(grep("Superior", classes)) > 1){
          choice <- sample(probes[grepl("Superior", classes)], 1)
        } else {
          choice <- probes[grepl("Superior", classes)]
        }
      }
      else if (any(grepl("Inferior", classes))){
        if(all(grepl("Inferior", classes))){
          choice <- sample(probes, 1)
        } else {
          choice <- sample(probes[!grepl("Inferior", classes)], 1)
        }
      }
      else {
        choice <- sample(probes, 1)
      }
      means <- object[choice,]
      retain <- choice
      dups <- probes[!probes==choice]
      return(list(means, retain, dups))
    }
    )
    means <- do.call("rbind", lapply(precchoice, function (x) x[[1]]))
    rownames(means) <- unlist(lapply(precchoice, function (x) x[[2]]))
    object[rownames(means),] <- means
    dups <- unlist(lapply(precchoice, function (x) x[[3]]))
    EPICv2manifest <- EPICv2manifest[!rownames(EPICv2manifest) %in% dups,]
    object <- object[!rownames(object) %in% dups,]
  }, random = {
    message("Selecting replicate probes at random...")
    randchoice <- lapply(posreps, function (x) {
      probes <- rownames(EPICv2manifest)[coords==x]
      choice <- sample(probes, 1)
      dups <- probes[!probes==choice]
      list(choice, dups)
    })
    remove <- unlist(lapply(randchoice, function (x) x[[2]]))
    EPICv2manifest <- EPICv2manifest[!rownames(EPICv2manifest) %in% remove,]
    object <- object[!rownames(object) %in% remove,]
    
  })
  object
}
