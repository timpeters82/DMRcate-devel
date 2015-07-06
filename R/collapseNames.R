## for tidying up gene names and associations
collapseNames <- function(x, glue=';', coll=',') {  ## another helper

  uunlist <- function(x) unique(unlist(x)) # helper for common idiom
  if(!is.character(x)) x <- as.character(x)
  if(length(x) > 1) {
    x <- paste(do.call(c, lapply(x, collapseNames)), collapse=glue)
  }
  if(any(grepl(glue, as.character(x)))) {
    x <- paste(uunlist(strsplit(as.character(x), glue)), collapse=coll)
  }
  if(grepl(coll, x)) {
    x <- paste(unique(strsplit(x, coll)[[1]]), collapse=coll)
  }
  return(x)

}
