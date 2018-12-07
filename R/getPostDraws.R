#' get posterior draws from a model object
#'
#' @param runObject an object returned by runJAGSmodel.
#' @return a list of model parameter posterior draws.
#' @examples
#' x <- getPostDraws(runObject)
#' @export
getPostDraws <- function(x) x$JAGSout$BUGSoutput$sims.list
