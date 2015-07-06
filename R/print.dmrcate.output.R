#' make it harder for me to shoot myself in the foot calling `some_DMR_results`
#'  
#' @param an object of class dmrcate.output
#'
#' @export
print.dmrcate.output <- function(object) { 
  cat("Object of class", class(object), "\n\n")
  cat("  Input data:", nrow(object$input), "features\n")
  cat("     Results:", nrow(object$results), "DMRs\n")
  cat("      Cutoff: p < ", object$cutoff, "\n")
}
