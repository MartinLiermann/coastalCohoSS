#' get table summarizing the posterior distributions for the model parameters
#'
#' @param runObject an object returned by runJAGSmodel.
#' @return a table summarizing the model parameter posterior distributions.
#' @examples
#' x <- getPosteriorSummary(runObject)
#' @export
getPosteriorSummary <- function(runObject){
  runObject$JAGSout$BUGSoutput$summary
}
