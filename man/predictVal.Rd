% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictVal.R
\name{predictVal}
\alias{predictVal}
\title{Generate model predictions based on the posterior}
\usage{
predictVal(runObj, predType = "smoltAtTrap", at = NULL)
}
\arguments{
\item{runObj}{an object created by runJAGSmodel}

\item{predType}{the value to be predicted}

\item{at}{the x locations at which predictions should be made}
}
\value{
a 3 dimensional array (sims, population, at value)
}
\description{
Generate model predictions based on the posterior
}
\examples{
predictVal(m7,"smoltAtDat",seq(0,10000,by=100))
}
