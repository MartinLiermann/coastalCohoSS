#' plots two variables against each other
#'
#' @param xx the x value (either a vector or a matrix where rows represent the MCMC sims).
#' @param yy the y value (see above).
#' @param xUn if TRUE a credible interval will be plotted for the x variable.
#' @param yUn if TRUE a credible interval will be plotted for the y variable.
#' @param xLab the x label.
#' @param yLab the y label.
#' @param xLim the x limits
#' @param yLim the y limits
#' @param main the main text
#' @param p the probability used to define the credible interval.
#' @param pch the point type passed to plot.
#' @return creates a plot.
#' @examples
#' plotXY(x,y)
#' @export
plotXY <- function(xx,yy,xUn=FALSE,yUn=FALSE,xLab="",yLab="",
                   xLim=c(NA,NA),yLim=c(NA,NA),main="",p=0.8,pch=16){

  minP <- function(x) if(is.null(dim(x))) min(x,na.rm=T) else min(apply(x,2,quantile,prob=0.25,na.rm=T))
  maxP <- function(x) if(is.null(dim(x))) max(x,na.rm=T) else max(apply(x,2,quantile,prob=0.75,na.rm=T))

  if(is.na(xLim[1])) xLim[1] <- minP(xx)
  if(is.na(xLim[2])) xLim[2] <- maxP(xx)
  if(is.na(yLim[1])) yLim[1] <- minP(yy)
  if(is.na(yLim[2])) yLim[2] <- maxP(yy)

  pCol <- rgb(0.1,0.9,0.1,0.5)
  pCex <- 1.5

  # check the dimension of xx and yy and calculate plotting values
  x2D <- !is.null(dim(xx))
  y2D <- !is.null(dim(yy))
  if(x2D) xVals <- apply(xx,2,median) else xVals <- xx
  if(y2D) yVals <- apply(yy,2,median) else yVals <- yy

  plot(xVals,yVals,bty="l",pch=pch,xlim=xLim,ylim=yLim,xlab=xLab,ylab=yLab,
       main=main,xaxs="r",yaxs="r",cex=pCex,col=pCol)

  # plot p% credible intervals for quantities if indicated and possible
  if(xUn & x2D){
    xQ <- apply(xx,2,quantile,prob=c((1-p)/2,(1+p)/2))
    segments(x0=xQ[1,],y0=yVals,x1=xQ[2,],y1=yVals,lwd=3,col=pCol)
  }
  if(yUn & y2D){
    yQ <- apply(yy,2,quantile,prob=c((1-p)/2,(1+p)/2))
    segments(y0=yQ[1,],x0=xVals,y1=yQ[2,],x1=xVals,lwd=3,col=pCol)
  }
  points(xVals,yVals,pch=pch,cex=pCex,col=pCol)
}


