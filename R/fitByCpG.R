## convenience function for comparison purposes & arbitrary marginal test stats
##
fitByCpG <- function(x, designOrSurv, what=c('linear','cox','logistic')) { 

  what <- match.arg(what)

  if (what == 'linear') {

    require(limma)
    xx <- getM(x)
    if (any(is.na(xx))) {
      require(impute)
      xx <- impute.knn(xx)$data
    }
    fit <- lmFit(xx, design)
    return(eBayes(fit)) ## should return the T-statistic instead

  } else if (what == 'cox') { 

    require(simulatorZ)

  } else { 
    stop(paste(what, 'is not supported as a response yet.'))
  }

} 
