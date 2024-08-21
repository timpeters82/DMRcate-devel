cpg.annotate <- function (datatype = c("array", "sequencing"), object, what = c("Beta", "M"), 
          arraytype = c("EPICv2", "EPICv1", "EPIC", "450K"), epicv2Remap = TRUE,
          epicv2Filter = c("mean", "sensitivity", "precision", "random"),
          analysis.type = c("differential", "variability", "ANOVA", "diffVar"), 
          design, contrasts = FALSE, cont.matrix = NULL, fdr = 0.05, coef, 
          varFitcoef = NULL,  topVarcoef = NULL, ...) 
{
  analysis.type <- match.arg(analysis.type)
  what <- match.arg(what)
  arraytype <- match.arg(arraytype)
  if (arraytype == "EPIC") {
    stop("Please specify either 'EPICv2' or 'EPICv1' for arraytype. EPICv2 probe IDs have 15 characters, e.g. cg00000029_TC21. EPICv1 probe IDs have 10 characters, e.g. cg00000029.")
  }
  if (datatype == "array") {
    stopifnot(class(object)[1] %in% c("matrix", "GenomicRatioSet"))
    if (is(object, "matrix")) {
      #Only retain cytosines
      keep <- grep("^cg|^ch", rownames(object))
      object <- object[keep,]
        if (arraytype == "450K") {
        grset <- makeGenomicRatioSetFromMatrix(mat = object, 
                                               array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", 
                                               what = what)
        
      }
      if (arraytype == "EPICv1") {
        grset <- makeGenomicRatioSetFromMatrix(mat = object, 
                                               array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19", 
                                               what = what)
        
      }
      if (arraytype == "EPICv2") {
        grset <- makeGenomicRatioSetFromMatrix(mat = object, 
                                               array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38", 
                                               what = what)
        anno <- getAnnotation(grset)
        message("EPICv2 specified. Loading manifest...")
        ah <- AnnotationHub()
        EPICv2manifest <- ah[["AH116484"]]
        anno <- cbind(anno, EPICv2manifest[rownames(anno), 73:80])
        object <- getM(grset)
      } 
      
    } else {
    grset <- object
    anno <- getAnnotation(grset)
    }
    
    if(arraytype!="EPICv2"){
      object <- getM(grset)
    } else {
      #Remapping
      if(epicv2Remap){
        if(any(anno$CH_WGBS_evidence=="Y")){
          #Remap those with offtarget
          torm <- sum(anno$CH_WGBS_evidence=="Y")
          message(paste0("Remapping ", torm, " cross-hybridising probes to their more likely offtarget..."))
          anno$chr[anno$CH_WGBS_evidence=="Y"] <- gsub(":.*", "", anno$Suggested_offtarget[anno$CH_WGBS_evidence=="Y"])
          anno$pos[anno$CH_WGBS_evidence=="Y"] <- as.integer(gsub(".*:", "", anno$Suggested_offtarget[anno$CH_WGBS_evidence=="Y"]))
          #Throw away those with CH but no remap
          anno <- anno[!(anno$CH_BLAT=="Y" & anno$CH_WGBS_evidence==""),]
          object <- object[rownames(anno),]
        }
      }
      
      #Check for replicates
      coords <- paste(anno$chr, anno$pos, sep=":")
      posreps <- table(coords)
      if (any(posreps > 1)){
        message(paste("Replicate probes that map to the same CpG site found. Filtering these using strategy:", epicv2Filter))
        if(any(nchar(rownames(object)) < 13)){
          stop("Error: rownames do not look like EPICv2 probes. This function will only work for EPICv2 data.")
        }
        if(any(!grepl("^cg|^ch", rownames(object)))){
          stop("Error: This function will only accept a matrix with rownames beginning with cg or ch. Please run your matrix through rmSNPandCH() first.")
        }
        posreps <- names(posreps)[posreps > 1]
        switch(epicv2Filter, mean={
          message("Averaging probes that map to the same CpG site...")
          outs <- lapply(posreps, function (x){
            ids <- coords==x
            means <- colMeans(object[ids,])
            retain <- rownames(anno)[ids][1]
            dups <- rownames(anno)[ids][-1]
            list(means, retain, dups)
          })
          means <- do.call("rbind", lapply(outs, function (x) x[[1]]))
          rownames(means) <- unlist(lapply(outs, function (x) x[[2]]))
          object[rownames(means),] <- means
          dups <- unlist(lapply(outs, function (x) x[[3]]))
          anno <- anno[!rownames(anno) %in% dups,]
          object <- object[!rownames(object) %in% dups,]
        }, sensitivity={
          message("Selecting probes that map to the same CpG site by sensitivity to methylation change...")
          senschoice <- lapply(posreps, function (x) {
            probes <- rownames(anno)[coords==x]
            classes <- anno[probes, "Rep_results_by_LOCATION"]
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
          anno <- anno[!rownames(anno) %in% remove,]
          object <- object[!rownames(object) %in% remove,]
        }, precision={
          message("Processing probes that map to the same CpG site for best precision...")
          precchoice <- lapply(posreps, function (x) {
            probes <- rownames(anno)[coords==x]
            classes <- anno[probes, "Rep_results_by_LOCATION"]
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
          anno <- anno[!rownames(anno) %in% dups,]
          object <- object[!rownames(object) %in% dups,]
        }, random = {
          message("Selecting replicate probes at random...")
          randchoice <- lapply(posreps, function (x) {
            probes <- rownames(anno)[coords==x]
            choice <- sample(probes, 1)
            dups <- probes[!probes==choice]
            list(choice, dups)
          })
          remove <- unlist(lapply(randchoice, function (x) x[[2]]))
          anno <- anno[!rownames(anno) %in% remove,]
          object <- object[!rownames(object) %in% remove,]
          
        })
      }
    }
    switch(analysis.type, differential = {
      stopifnot(is.matrix(design))
      if (!contrasts) {
        stopifnot(colnames(design)[1] == "(Intercept)")
      } else {
        stopifnot(!is.null(cont.matrix))
      }
      fit <- lmFit(object, design, ...)
      if (contrasts) {
        stopifnot(coef %in% colnames(cont.matrix))
        fit <- contrasts.fit(fit, cont.matrix)
      }
      fit <- eBayes(fit)
      tt <- topTable(fit, coef = coef, number = nrow(object))
      nsig <- sum(tt$adj.P.Val < fdr, na.rm = T)
      if (nsig == 0) {
        message("Your contrast returned no individually significant probes. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DMRs, but be warned there is an increased risk of Type I errors.")
      }
      if (nsig > 0 & nsig <= 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes; a small but real effect. Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases the risk of Type I errors."))
      }
      if (nsig > 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes. We recommend the default setting of pcutoff in dmrcate()."))
      }
      betafit <- lmFit(ilogit2(object), design, ...)
      if (contrasts) {
        betafit <- contrasts.fit(betafit, cont.matrix)
      }
      betafit <- eBayes(betafit)
      betatt <- topTable(betafit, coef = coef, number = nrow(object))
      m <- match(rownames(tt), rownames(betatt))
      tt$diff <- betatt$logFC[m]
      m <- match(rownames(object), rownames(tt))
      tt <- tt[m, ]
      stat <- tt$t
      rawpval <- tt$P.Value
      
      annotated <- GRanges(paste(anno$chr, anno$pos, sep=":"),
                             stat = stat, diff = tt$diff, rawpval = rawpval,
                             ind.fdr = tt$adj.P.Val, 
                             is.sig = tt$adj.P.Val < fdr)
        
      names(annotated) <- paste(annotated)
      annotated <- sort(sortSeqlevels(annotated))
      
    }, variability = {
      wholevar <- var(object)
      weights <- apply(object, 1, var)
      weights <- weights/mean(weights)
      annotated <- GRanges(as.character(anno$chr), IRanges(anno$pos, 
                                                             anno$pos), stat = weights, diff = rep(0, nrow(object)), 
                           rawpval = rep(0, nrow(object)), ind.fdr = rep(0, nrow(object)), is.sig = weights > 
                             quantile(weights, 0.95))
      names(annotated) <- rownames(object)
    }, ANOVA = {
      message("You are annotating in ANOVA mode: consider making the value of fdr quite small, e.g. 0.001")
      stopifnot(is.matrix(design))
      fit <- lmFit(object, design, ...)
      fit <- eBayes(fit)
      sqrtFs <- sqrt(fit$F)
      sqrtfdrs <- p.adjust(fit$F.p.value, method = "BH")
      rawpval <- fit$F.p.value
      nsig <- sum(sqrtfdrs < fdr)
      if (nsig == 0) {
        message("Your design returned no individually significant probes for ANOVA. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DMRs, but be warned there is an increased risk of Type I errors.")
      }
      if (nsig > 0 & nsig <= 100) {
        message(paste("Your design returned", nsig, 
                      "individually significant probes for ANOVA; a small but real effect. Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases the risk of Type I errors."))
      }
      if (nsig > 100) {
        message(paste("Your design returned", nsig, 
                      "individually significant probes for ANOVA. We recommend the default setting of pcutoff in dmrcate(). Large numbers (e.g. > 100000) may warrant a smaller value of the argument passed to fdr"))
      }
      anno <- getAnnotation(grset)
      stat <- sqrtFs
      annotated <- GRanges(as.character(anno$chr), IRanges(anno$pos, 
                                                           anno$pos), stat = stat, diff = 0, rawpval = rawpval, ind.fdr = sqrtfdrs, 
                           is.sig = sqrtfdrs < fdr)
      names(annotated) <- rownames(object)
    }, diffVar = {
      stopifnot(is.matrix(design))
      if (!contrasts) {
        stopifnot(colnames(design)[1] == "(Intercept)")
      } else {
        stopifnot(!is.null(cont.matrix))
      }
      fitvar <- varFit(object, design = design, coef = varFitcoef)
      if (contrasts) {
        fitvar <- contrasts.varFit(fitvar, cont.matrix)
      }
      tt <- topVar(fitvar, coef = topVarcoef, number = nrow(object))
      nsig <- sum(tt$Adj.P.Value < fdr)
      if (nsig == 0) {
        message("Your contrast returned no individually significant probes. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DVMRs, but be warned there is an increased risk of Type I errors.")
      }
      if (nsig > 0 & nsig <= 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes; a small but real effect. Consider manually setting the value of pcutoff to return more DVMRs, but be warned that doing this increases the risk of Type I errors."))
      }
      if (nsig > 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes. We recommend the default setting of pcutoff in dmrcate()."))
      }
      m <- match(rownames(object), rownames(tt))
      tt <- tt[m, ]
      anno <- getAnnotation(grset)
      stat <- tt$t
      rawpval <- tt$P.Value
      annotated <- GRanges(as.character(anno$chr), IRanges(anno$pos, 
                                                           anno$pos), stat = stat, diff = 0, 
                           rawpval=rawpval, ind.fdr = tt$Adj.P.Value, 
                           is.sig = tt$Adj.P.Value < fdr)
      names(annotated) <- rownames(tt)
    })
    annotated <- sort(annotated)
    return(new("CpGannotated", ranges = annotated))
  }
  if (datatype == "sequencing") {
    stop("Sequencing mode is deprecated for cpg.annotate(). Please use sequencing.annotate().")
  }
  else {
    message("Error: datatype must be one of 'array' or 'sequencing'")
  }
}
