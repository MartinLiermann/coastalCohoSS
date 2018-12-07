#' a function to retrieve the Rhat values from the model object
#'
#' @param runObject an object returned by runJAGSmodel.
#' @return the Rhat values for the model parameters.
#' @examples
#' rr <- getRhat(runObject)
#' @export
getRhat <- function(runObject){
  runObject$JAGSout$BUGSoutput$summary[,"Rhat"]
}
