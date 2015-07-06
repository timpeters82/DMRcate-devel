##
## Convenience and filtering functions for fitting mixture models 
## (such as the shifted reflected Gamma originally implemented by Pan Du, or 
## Raftery's mclust) to DMRs across samples rather than sites within a sample.
##
## First step: Fitting a mixture model via E-M is not terribly fast,
## so it helps to get rid of sites/DMRs where there is very little chance
## of finding a useful mixture model fit at all. We'll go back and re-fit
## the entirety of DMRs across the entirety of a study later, after tidying.
##
mixFilterDMRs <- function(grset, maxMinM=-2, minMaxM=2, minRange=6) {  # {{{

  M <- prepM(grset)
  maxM <- rowMaxs(M, na.rm=T)
  minM <- rowMins(M, na.rm=T)
  rngM <- maxM - minM

  keep <- which( maxM >= minMaxM & minM <= maxMinM & rngM >= minRange)
  message('Keeping ', length(keep), ' out of ', nrow(grset), ' DMRs...')

  trimmed <- grset[ keep, ]
  metadata(trimmed)$mixFiltered <- TRUE
  return(trimmed)

} # }}}


## Check before starting to fit mixtures that the fits are likely to succeed
##
isMixFiltered <- function(grset) { # {{{
  ('mixFiltered' %in% names(metadata(grset))) && (metadata(grset)$mixFiltered)
} # }}}


## Next, fit a mixture model to the filtered subset of global DMRs
## (the reason we default to a Gaussian can be seen using plotGammaFit)
##
fitDmrMixture <- function(grset, how=c('gaussian', 'gamma')) { # {{{

  if (!isMixFiltered(grset)) grset <- mixFilterDMRs(grset)
  bigVec <- as.vector(prepM(grset))
  how <- match.arg(how)
  switch(how, 
         gaussian=suppressWarnings(Mclust(bigVec, G=1:3)), 
         gamma=gammaFitEM(bigVec))

} # }}}


## Given a mixture fit, generate a DMR barcode "caller" for a dataset.  
##
## I have found that it's a good idea to "tighten" calls of partial methylation,
## which for my purposes means "possibly imprinted or allelically methylated". 
## You my not want to do this, in which case set maxProbPartial=0.5 or higher.
##
## Right now this only works for Gaussian mixture models (because Gamma didn't
## do a good job on DMR sets that contain substantial partial methylation).
##
generateCaller <- function(mixfit, maxProbPartial=.05) { # {{{
  
  if (class(mixfit) != 'Mclust') 
    stop('Only Gaussian mixture models are supported at the moment.')

  values <- seq(from=0, to=1, length.out=mixfit$G)
  mixparams <- data.frame(mean=mixfit$parameters$mean,
                          sd=sqrt(mixfit$parameters$variance$sigmasq))
  rownames(mixparams) <- c('unmeth','partial','meth')

  dunmeth <- function(x) with(mixparams['unmeth',], dnorm(x, mean=mean, sd=sd))
  dmeth <- function(x) with(mixparams['meth',], dnorm(x, mean=mean, sd=sd))
  if (mixfit$G > 2)
    dpartial <- function(x) with(mixparams['partial',],dnorm(x,mean=mean,sd=sd))

  barcodeDMR <- function(dmr, withZ=F) {
    if (class(dmr) %in% c('matrix','data.frame','DataFrame')) {
      res <- t(apply(dmr, 1, barcodeDMR))
      colnames(res) <- colnames(dmr)
    } else { 
      if (mixfit$G == 2) {
        Z <- cbind(dunmeth(dmr), dmeth(dmr))
      } else { 
        Z <- cbind(dunmeth(dmr), dpartial(dmr), dmeth(dmr))
      }
      tighten <- function(zz) {
        if (mixfit$G == 3 && which.max(zz) == 2) {
          if (zz[1] > maxProbPartial) return(1)
          else if (zz[3] > maxProbPartial) return(3)
          else return(2)
        } else return(which.max(zz))
      }
      res <- values[ apply(Z, 1, tighten) ]
      if (withZ == TRUE) attr(res, 'Z') <- Z
    }
    return(res)
  }
  attr(barcodeDMR, 'mixparams') <- mixparams
  return(barcodeDMR)

} # }}}


## Finally, fit the resulting caller to the entire dataset of DMRs. 
## 
callDMRs <- function(grset, caller=NULL) { # {{{

  if (is.null(caller)) caller <- generateCaller(fitDmrMixture(grset))
  res <- t(apply(prepM(grset), 1, caller))
  colnames(res) <- colnames(grset)
  attr(res, 'mixparams') <- attr(caller, 'mixparams')
  attr(res, 'caller') <- caller
  return(res)

} # }}}
