loctest <- function (fhat1, fhat2) 
    ## Adapted from ks:::kde.local.test.1d. Takes in 2 density
    ## estimates fhat1 and fhat2 and provides a local signficance test
    ## between the 2 estimates at each of fhat1$eval.points ==
    ## fhat2$eval.points.  fhat1$eval.points will always be identical
    ## to fhat2$eval.points given that dmrcate() supplies the same
    ## argument for both in kde()
{
    n1 <- length(fhat1$x)
    n2 <- length(fhat2$x)
    d <- 1
    RK <- (4 * pi)^(-d/2)
    h1 <- fhat1$h
    h2 <- fhat2$h
    h2D2fhat <- 0
    fhat.diff <- fhat1$estimate - fhat2$estimate - 1/2 * h2D2fhat
    var.fhat.diff <- ((n1 * h1)^(-1) * fhat1$estimate +
                      (n2 * h2)^(-1) * fhat2$estimate) * RK
    X2 <- fhat.diff^2/var.fhat.diff
    pvalue <- pchisq(X2, 1, lower.tail=FALSE)
                                        # Changed from 'pvalue <- 1 -
                                        # pchisq(X2, 1)' to get values
                                        # < 2.2e-16
    result <- list(fhat1 = fhat1, fhat2 = fhat2, chisq = X2, 
                   pvalue = pvalue)
    class(result) <- "kde.loctest"
    return(result)
}

annotate <- function(object,
    annotation=c(array="IlluminaHumanMethylation450k",
      annotation="ilmn12.hg19"),
    analysis.type=c("differential", "variability", "hypermethylation"),
    design, coef, diff.metric=c("FC", "t"), paired=FALSE, pcutoff=NULL,
    betacutoff=NULL, quantcut=0.9, cut=0.5,...)
{
    stopifnot(is.matrix(object))
    analysis.type <- match.arg(analysis.type)
    switch(analysis.type, "differential"={
        stopifnot(is.matrix(design))
        fit <- lmFit(object, design, ...)
        if (!paired) {
            fit <- contrasts.fit(fit, c(1, -1))
        }
        fit <- eBayes(fit)
        tt <- topTable(fit, coef=coef, number=nrow(object))
        if (!is.null(pcutoff) & !is.null(betacutoff)) {
            txt <- "Cannot specify both a p-value and beta threshold; 
                    one or neither needed"
            stop(paste(strwrap(txt, exdent=2), collapse="\n"))
        }
        if (!is.null(pcutoff)) {
            stopifnot(0 <= pcutoff & pcutoff <= 1)
            tt <- tt[tt$adj.P.Val <= pcutoff,]
        }
        if (!is.null(betacutoff)) {
            stopifnot(0 < betacutoff & betacutoff < 1)
            betafit <- lmFit(ilogit2(object), design, ...)
            if (!paired) {
                betafit <- contrasts.fit(betafit, c(1, -1))
            }
            betafit <- eBayes(betafit)
            betatt <- topTable(betafit, coef=coef, number=nrow(object))
            probesubs <- betatt[abs(betatt$logFC) > betacutoff,]
            tt <- tt[rownames(tt) %in% rownames(probesubs),]
        }
        m <- match(rownames(tt), rownames(object))
        object <- object[m,]
        RSobject <- RatioSet(object, annotation=annotation)
        RSanno <- getAnnotation(RSobject)
        weights <- switch(diff.metric, "FC"={
            tt$logFC
        }, "t"={
            tt$t
        })
        annotated <- data.frame(ID=rownames(object), weights=weights,
            CHR=RSanno$chr, pos=RSanno$pos, gene=RSanno$UCSC_RefGene_Name,
            group=RSanno$UCSC_RefGene_Group)
    }, "variability"={
        stopifnot(0 <= quantcut & quantcut <= 1)
        RSobject <- RatioSet(object, annotation=annotation)
        RSanno <- getAnnotation(RSobject)
        annotated <- data.frame(ID=rownames(object),
            weights=apply(object, 1, var), CHR=RSanno$chr, pos=RSanno$pos,
            gene=RSanno$UCSC_RefGene_Name, group=RSanno$UCSC_RefGene_Group)
        if (quantcut > 0) {
            annotated <- annotated[annotated$weights >
                quantile(annotated$weights, quantcut),]
        }
    }, "hypermethylation"={
        stopifnot(0 <= cut & cut <= 1)
        RSobject <- RatioSet(object, annotation=annotation)
        RSanno <- getAnnotation(RSobject)
        annotated <- data.frame(ID=rownames(object),
            weights=apply(object, 1, mean), CHR=RSanno$chr, pos=RSanno$pos,
            gene=RSanno$UCSC_RefGene_Name, group=RSanno$UCSC_RefGene_Group)
        annotated <- annotated[annotated$weights > logit2(cut),]
    })
    annotated <- annotated[order(annotated$CHR, annotated$pos),]
    class(annotated) <- "annot"
    return(annotated)
}

rmSNPandCH <- function(object, dist=2, mafcut=0.05, and=TRUE, rmcrosshyb=TRUE)
{
    stopifnot(is.matrix(object))
    dist <- as.integer(dist)
    stopifnot(0 <= mafcut & mafcut <= 1)

    env <- new.env(parent=emptyenv())
    data(dmrcatedata, envir=env)
    len0 <- vapply(env$illuminaSNPs$Distance, length, integer(1),
                   USE.NAMES=FALSE)
    len1 <- vapply(env$illuminaSNPs$MinorAlleleFrequency, length,
                   integer(1), USE.NAMES=FALSE)
    keep <- len0 == len1
    len0 <- len0[keep]
    env$illuminaSNPs <- env$illuminaSNPs[keep,]

    dist0 <- as.integer(unlist(env$illuminaSNPs$Distance, use.names=FALSE))
    distrange <- range(dist0)
    stopifnot(dist >= min(distrange) && dist <= max(distrange))

    test0 <- (dist0 >= -1) & (dist0 <= dist)
    test1 <-
        unlist(env$illuminaSNPs$MinorAlleleFrequency, use.names=FALSE) > mafcut
    test <- if (and) (test0 & test1) else (test0 | test1)
    ## 'any' by group
    ntrue <- cumsum(test)[cumsum(len0)]
    badidxs <- ntrue - c(0, head(ntrue, -1)) != 0

    badprobes <- rownames(env$illuminaSNPs)[badidxs]
    if(rmcrosshyb)
        badprobes <- union(badprobes, as.character(env$crosshyb))

    object[!(rownames(object) %in% badprobes),]
}

dmrcate <- function(object, bw=1000, p.adjust.method="BH", pcutoff=0.05,
                    consec=FALSE, consecbw=10)
{
    stopifnot(is(object, "annot"))
    stopifnot(bw >= 1)
    stopifnot(0 <= pcutoff & pcutoff <= 1)
    object <- data.frame(ID=object$ID, weights=object$weights, CHR=object$CHR,
                         pos=object$pos, gene=object$gene, group=object$group)
    chromosome <- unique(c(as.character(object$CHR)))
    names(chromosome) <- unique(c(as.character(object$CHR)))
    if(consec & is.null(consecbw)){
        stop("Consecutive probe bandwidth must be specified")
    }
    for (i in chromosome)
    {
        message("Fitting ", chromosome[i], "...", sep='')
        chrIndex <- object$CHR %in% chromosome[i]
        if (consec) {
            stopifnot(consecbw >= 1)
            posIndex <- seq_along(object$pos[chrIndex])
            densest <- suppressWarnings({
                kde(posIndex, h=consecbw, eval.points=posIndex,
                    w=abs(object$weights[chrIndex]))
            })
            unweighted <- kde(posIndex, h=consecbw, eval.points=posIndex)
            test <- loctest(fhat1=densest, fhat2=unweighted)
        } else {
            pos <- object$pos[chrIndex]
            densest <- suppressWarnings({
                kde(pos, h=bw, eval.points=pos, w=abs(object$weights[chrIndex]))
            })
            unweighted <- kde(pos, h=bw, eval.points=pos)
            test <- loctest(fhat1=densest, fhat2=unweighted)
        }
        
        object$raw[chrIndex] <- test$pvalue
        
    }
    object$fdr <- p.adjust(object$raw, method=p.adjust.method)
    sigprobes <- object[object$fdr <= pcutoff,, drop=FALSE]
    if (nrow(sigprobes)==0) {
        txt <- "No signficant regions found. Try increasing the value of
                'pcutoff' in 'dmrcate' and/or 'annotate'."
        stop(paste(strwrap(txt, exdent=2), collapse="\n"))
    }
    probeIDs <- pvals <- genes <- group <- hg19coord <- list()
    message("Demarcating regions...")

    regions <- data.frame(ID=probeIDs, pvals=pvals, genesassoc=genes,
        group=group, hg19coord=hg19coord)

    ## coerce types
    sigprobes$ID <- as.character(sigprobes$ID)
    sigprobes$fdr <- as.numeric(sigprobes$fdr)
    sigprobes$gene <- as.character(sigprobes$gene)
    sigprobes$group <- as.character(sigprobes$group)
    for(i in unique(sigprobes$CHR)) {
        this.chr <- sigprobes[sigprobes$CHR == i,]
        this.chr$diffs <- c(Inf, diff(this.chr$pos))
        regstarts <- c(which(this.chr$diffs > bw), nrow(this.chr) + 1)
        nRegstarts1 <- length(regstarts) - 1
        ## pre-allocate space for results
        probeIDs <- pvals <- genes <- group <- hg19coord <-
            vector("list", nRegstarts1)
        for (j in seq_len(nRegstarts1)) {
            regIndex <- regstarts[j]:(regstarts[j+1]-1)
            probeIDs[j] <- list(this.chr$ID[regIndex])
            pvals[j] <- list(this.chr$fdr[regIndex])
            genes[j] <- list(this.chr$gene[regIndex])
            group[j] <- list(this.chr$group[regIndex])
            hg19coord[j] <-
                list(paste(chromosome[i], ":", this.chr$pos[regstarts[j]], "-",
                           this.chr$pos[regstarts[j+1]-1], sep=''))
        }
        regions <- rbind(regions, data.frame(
            probeIDs=as.data.frame(as.matrix(probeIDs)),
            pvals=as.data.frame(as.matrix(pvals)),
            genes=as.data.frame(as.matrix(genes)),
            group=as.data.frame(as.matrix(group)),
            hg19coord=as.data.frame(as.matrix(hg19coord))))
    }
    ## Assigning relevant annotation and p-values to regions
    uunlist <- function(x)              # helper for common idiom
        unique(unlist(x))
    if (nrow(regions) > 1) {
        results <- data.frame(
            gene_assoc=as.data.frame(as.matrix(sapply(regions[,3], uunlist))),
            group=as.data.frame(as.matrix(sapply(regions[,4], uunlist))),
            hg19coord=as.data.frame(as.matrix(sapply(regions[,5], uunlist))),
            no.probes=unlist(lapply(regions[,1], length)),
            minpval=unlist(lapply(regions[,2], min)),
            meanpval=unlist(lapply(regions[,2], mean)))
    } else if (nrow(regions) ==1) {
        results <- data.frame(
            gene_assoc=as.data.frame(as.matrix(lapply(regions[,3], uunlist))),
            group=as.data.frame(as.matrix(lapply(regions[,4], uunlist))),
            hg19coord=as.data.frame(as.matrix(lapply(regions[,5], uunlist))),
            no.probes=unlist(lapply(regions[,1], length)),
            minpval=unlist(lapply(regions[,2], min)),
            meanpval=unlist(lapply(regions[,2], mean)))
    }
    colnames(results) <- c("gene_assoc", "group", "hg19coord",
                           "no.probes", "minpval", "meanpval")
    if(sum(grepl(';', unlist(results$gene_assoc))) > 0) {
        results$gene_assoc <- lapply(results$gene_assoc, function (x) {
            paste(unique(unlist(strsplit(unlist(x), ";"))), collapse=',')
        })
    }  
    results$gene_assoc <- unlist(results$gene_assoc)
    if(sum(grepl(';', unlist(results$group))) > 0) {
        results$group <- lapply(results$group, function (x) {
            paste(unique(unlist(strsplit(unlist(x), ";"))), collapse=',')
        })
    }
    results$group <- unlist(results$group)
    results <- results[order(results$minpval, -results$no.probes),, drop=FALSE]
    message("Done!")
    output <- NULL
    output$input <- object
    output$results <- results
    output$cutoff <- pcutoff
    return(output)
}

makeBedgraphs <- function(dmrcoutput, betas,
      annotation=c(array="IlluminaHumanMethylation450k",
        annotation="ilmn12.hg19"),
      samps=NULL)
{
    stopifnot(attributes(dmrcoutput)$names == c("input", "results", "cutoff"))
    stopifnot(is.matrix(betas))
    stopifnot(all(samps %in% 1:ncol(betas)))
    if (is.null(samps)) {
        samps <- 1:ncol(betas)
    }
    sigprobes <- dmrcoutput$input[dmrcoutput$input$fdr < dmrcoutput$cutoff,]
    betas <- betas[as.character(sigprobes$ID),]
    RSobject <- RatioSet(betas, annotation=annotation)
    RSanno <- getAnnotation(RSobject)
    probes <- matrix(nrow=nrow(betas), ncol=3)
    ## Below is making sure the bedgraph takes into account the
    ## strandedness of the DNA that the probe hybridises to
    probes[RSanno$strand == "-",] <-
        cbind(as.character(sigprobes$CHR)[RSanno$strand == "-"],
              as.integer(as.character(sigprobes$pos))[RSanno$strand == "-"],
              as.integer(as.character(sigprobes$pos))[RSanno$strand == "-"] + 1)
    probes[RSanno$strand == "+",] <-
        cbind(as.character(sigprobes$CHR)[RSanno$strand == "+"],
              as.integer(as.character(sigprobes$pos))[RSanno$strand == "+"] - 1,
              as.integer(as.character(sigprobes$pos))[RSanno$strand == "+"])
    for(i in samps) {
        write.table(cbind(probes, betas[,i]),
            paste(colnames(betas)[i], ".bedGraph", sep=''), sep='\t',
            row.names=FALSE, col.names=FALSE, quote=FALSE)
    } 
    message(length(samps), " bedGraphs written to ", getwd(), ".", sep='')
}

DMR.plot <- function(dmrcoutput, dmr, betas, phen.col,
    annotation=c(array="IlluminaHumanMethylation450k",
      annotation="ilmn12.hg19"),
    samps=NULL, toscale=FALSE, ...)
{
    stopifnot(attributes(dmrcoutput)$names == c("input", "results", "cutoff"))
    stopifnot(is.matrix(betas))
    stopifnot((length(dmr) == 1) && (dmr %in% 1:nrow(dmrcoutput$results)))
    stopifnot(all(samps %in% 1:ncol(betas)))
    stopifnot(ncol(betas) == length(phen.col))
    if(dmrcoutput$results$no.probes[dmr] < 2){
        stop("Region must have 2 or more CpGs")
    }
    coords <- dmrcoutput$results$hg19coord[dmr]
    chr <- sub(":.*", "", coords)
    bookends <- sub(".*:", "", coords)
    startcpg <- as.integer(sub("-.*", "", bookends))
    stopcpg <- as.integer(sub(".*-", "", bookends))
    RSobject <- RatioSet(betas, annotation=annotation)
    RSanno <- getAnnotation(RSobject)
    cpgs <- rownames(RSanno)[RSanno$chr %in% chr &
                             RSanno$pos >= startcpg & RSanno$pos <= stopcpg]
    cpgs <- cpgs[order(RSanno[cpgs,"pos"])]
    if(is.null(samps)){samps <- 1:ncol(betas)}
    betas <- betas[as.character(cpgs), samps]
    m <- match(as.character(cpgs), rownames(RSanno))
    clusta <- data.frame(
        gene=RSanno$UCSC_RefGene_Name[m],
        group=RSanno$UCSC_RefGene_Group[m], pos=RSanno$pos[m])
    if(!toscale) {
        plot(1:nrow(clusta), , type='n', xlab=paste(chr, "consecutive probes"),
             ylab="beta values", ylim=c(0,1.15), ...)
        for(i in 1:nrow(clusta)) {
            points(rep(i, ncol(betas)), betas[i,], col=phen.col, ...)
        } 
        empty <- rep(-1, nrow(clusta))
        empty[grep("5'UTR", clusta$group)] <- 1
        points(1:length(empty), empty, adj=0.6, col='blue4', pch=15)
        empty <- rep(-1, nrow(clusta))
        empty[grep("3'UTR", clusta$group)] <- 1.01
        points(1:length(empty), empty, adj=0.8, col='cyan', pch=15)
        empty <- rep(-1, nrow(clusta))
        empty[grep("Body", clusta$group)] <- 1.02
        points(1:length(empty), empty, col='firebrick1', pch=15)
        empty <- rep(-1, nrow(clusta))
        empty[grep("TSS1500", clusta$group)] <- 1.03
        points(1:length(empty), empty, adj=0.2, col='chartreuse', pch=15)
        empty <- rep(-1, nrow(clusta))
        empty[grep("TSS200", clusta$group)] <- 1.04
        points(1:length(empty), empty, adj=0.4, col='forestgreen', pch=15)
        empty <- rep(-1, nrow(clusta))
        empty[grep("1stExon", clusta$group)] <- 1.05
        points(1:length(empty), empty, adj=1, col='magenta', pch=15)
        empty <- rep('', nrow(clusta))
        geneidxs <- cumsum(rle(as.character(clusta$gene))$lengths)
        geneidxs <- ceiling(c(geneidxs[1]/2, geneidxs[2:length(geneidxs)] -
                              diff(geneidxs)/2))
        empty[geneidxs] <- as.character(clusta$gene[geneidxs])
        genevector <- empty
        text(1:nrow(clusta), 1.03, genevector, pos=3, offset=1)
        abline(0.5, 0, col='blue')
    }
    if (toscale) {
        plot(clusta$pos, 1:nrow(clusta), type='n',
             xlab=paste(chr, "hg19 coords"), ylab="beta values",
             ylim=c(0,1.15), ...)
        for(i in 1:nrow(clusta)){
            points(rep(clusta$pos[i], ncol(betas)), betas[i,], col=phen.col, ...)
        } 
        empty <- rep(-1, nrow(clusta))
        empty[grep("5'UTR", clusta$group)] <- 1
        points(clusta$pos, empty, adj=0.6, col='blue4', pch=15)
        empty <- rep(-1, nrow(clusta))
        empty[grep("3'UTR", clusta$group)] <- 1.01
        points(clusta$pos, empty, adj=0.8, col='cyan', pch=15)
        empty <- rep(-1, nrow(clusta))
        empty[grep("Body", clusta$group)] <- 1.02
        points(clusta$pos, empty, col='firebrick1', pch=15)
        empty <- rep(-1, nrow(clusta))
        empty[grep("TSS1500", clusta$group)] <- 1.03
        points(clusta$pos, empty, adj=0.2, col='chartreuse', pch=15)
        empty <- rep(-1, nrow(clusta))
        empty[grep("TSS200", clusta$group)] <- 1.04
        points(clusta$pos, empty, adj=0.4, col='forestgreen', pch=15)
        empty <- rep(-1, nrow(clusta))
        empty[grep("1stExon", clusta$group)] <- 1.05
        points(clusta$pos, empty, adj=1, col='magenta', pch=15)
        empty <- rep('', nrow(clusta))
        geneidxs <- cumsum(rle(as.character(clusta$gene))$lengths)
        geneidxs <- ceiling(c(geneidxs[1]/2, geneidxs[2:length(geneidxs)] -
                              diff(geneidxs)/2))
        empty[geneidxs] <- as.character(clusta$gene[geneidxs])
        genevector <- empty
        text(clusta$pos, 1.03, genevector, pos=3, offset=1)
        abline(0.5, 0, col='blue') 
    }
}
