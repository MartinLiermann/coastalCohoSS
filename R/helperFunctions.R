#' logit function
#'
#' @param x a real value
#' @return the logit of x
#' @export
logit <- function(x) log(x/(1-x))

#' inverse logit function
#'
#' @param x a real value
#' @return the inverse logit of x
#' @export
invLogit <- function(x) 1/(1+exp(-x))

#' fills gaps in a series
#'
#' @param vals the values which you would like to fill.
#' @param groups the variable which defines groups.
#' @param returnEstType if TRUE, the function returns the method used to fill.
#' @return a vector of values with gaps filled.
#' @export
fillGaps <- function(vals,groups,returnEstType=FALSE){
  estM <- rep(NA,length(vals))
  seVals <- tapply(vals,list(groups),mean,na.rm=T)
  seProx <- quantile(seVals,probs=c(0.75),na.rm=T)
  for(gg in unique(groups)){
    if(all(is.na(vals[groups==gg]))){
      vals[groups==gg] <- seProx
      estM[groups==gg] <- "otherPops"
    } else if(any(is.na(vals[groups==gg]))){
      ind <- is.na(vals[groups==gg])
      vals[groups==gg][ind] <- seVals[gg]
      estM[groups==gg][ind] <- "otherYrs"
    }
  }
  if(returnEstType){
    retVal <- list(vals=vals,estM=estM)
  }else{
    retVal <- vals
  }
}

#' fills gaps in a series using neighbors
#'
#' @param dat the JAGS data object.
#' @param var the variable to fill.
#' @param fKern the kernal used to fill.
#' @return a vector of values with gaps filled.
#' @export
fillGapsKern <- function(dat,param,fKern=c(rep(0.01,1000),1:3,0,3:1,rep(0.01,1000))){
  fKernInd <- (1:length(fKern)) - ceiling(length(fKern)/2)
  tmp <- dat[[param]]
  for(pop in 1:dat$Npops){
    allInds <- which(dat$stock==pop)
    naInds <- which(is.na(dat[[param]]) & dat$stock==pop)
    yrs <- dat$year[dat$stock==pop]-min(dat$year[dat$stock==pop])+1
    for(ind in naInds){
      useInds <- which((ind+fKernInd) %in% allInds)
      useInds <- useInds[!is.na(dat[[param]][ind + fKernInd[useInds]])]
      tmp[ind] <- sum(fKern[useInds]*dat[[param]][ind + fKernInd[useInds]])/
                  sum(fKern[useInds])
    }
  }
  tmp
}


#' Expands data to include all years and observations.
#'
#' @param data the data you would like to fill.
#' @param index the index indicating which of the total values were observed.
#' @param N the total number of length of the vector to be created.
#' @return a vector NAs of length N, with values at indices, index, replaced with data.
#' @export
expandData <- function(data,index,N){
  exDat <- rep(NA,N)
  exDat[index] <- data
  exDat
}

#' Expands data to include all years and observations for all data in bdat.
#'
#' @param bdat a list with the data you would like to fill.
#' @return a dataframe with all of the expanded data.
#' @export
expandAllData <- function(jdat){
  dataToExpand <- list(
    c("escapementObs","escInd"),
    c("escapementFemObs","escFemInd"),
    c("smoltObs","smoltInd","pHOS")
  )
  N <- jdat$N
  exData <- list(stock=jdat$stock,year=jdat$year)
  nn <- 3
  for(dp in dataToExpand){
    exData <- c(exData,list(expandData(jdat[dp[1]][[1]],jdat[dp[2]][[1]],N)))
    names(exData)[nn] <- dp[1]
    nn <- nn+1
  }
  for(nam in names(jdat)) if(!(nam %in% names(exData))) exData <- c(exData,jdat[nam])
  exData
}

#' Plot posteriorsDists.
#'
#' @param parSims a array with samples from the posterior distribution.
#' @param nSims the number of simulations in each chain
#' @return a plot of the posterior.
#' @export
plotPost <- function(parSims,parName="par",mainText="",nSims=1000){
  dd <- dim(parSims)
  nChains <- dd[1]/nSims
  n <- dd[2]
  qVals <- apply(parSims,2,quantile,prob=c(0.1,0.5,0.9))
  yRange <- range(qVals)
  plot(1,1,xlim=c(0,n),ylim=yRange,xlab="Index",ylab=parName,main=mainText,bty="l",type="n",xaxt="n")
  axis(side=1,at=(1:n)-0.5,labels=1:n)
  for(i in 1:nChains){
    qVals <- apply(parSims[1:nSims + (i-1)*nSims,],2,quantile,prob=c(0.1,0.5,0.9))
    segments(x0=(1:n-0.9+0.2*i),x1=(1:n-0.9+0.2*i),y0=qVals[1,],y1=qVals[3,],lwd=2,col=i)
    points((1:n-0.9+0.2*i),qVals[2,],pch=16,col=i)
  }
}

#' Plot predicted vs observed escapement.
#'
#' @param runObj an object returned by the runJAGSmodel function.
#' @return a plot of predicted vs observed escapement.
#' @export
plotEsc <- function(runObj,p=0.80,newWindow=FALSE,logScale=FALSE){
  qP <- (1-p)/2
  x <- getPostDraws(runObj)
  bdat <- runObj$dat$jagsDat
  odat <- runObj$dat$otherDat
  popNames <- aggregate(names(bdat$stock),list(bdat$stock),unique)$x
  if(newWindow) windows(8,8)
  par(mfrow=c(bdat$Npops,1),oma=c(3,3,1,1),mar=c(2,2,3,1))
  qEsc <- apply(x$wildEscapementAge4to5,2,quantile,prob=c(qP,0.5,1-qP))
  eInd <- bdat$escInd
  yrs <- odat$year
  esc <- rep(NA,bdat$N)
  esc[eInd] <- bdat$wildEscapementAge4to5obs
  yrRange <- range(odat$year)
  yLab <- "Escapement (4&5)"
  if(logScale){
    qEsc <- log(qEsc)
    esc <- log(esc)
    yLab <- "log Escapement (4&5)"
  }
  for(i in 1:bdat$Npops){
    pInd <- bdat$stock==i
    yRange <- c(0,max(qEsc[,pInd]))
    plot(yrs[pInd],esc[pInd],bty="l",xlab="", ylab="",pch=16,main=popNames[i],xlim=yrRange,ylim=yRange)
    lines(yrs[pInd],qEsc[2,pInd])
    lines(yrs[pInd],qEsc[1,pInd],lty=3)
    lines(yrs[pInd],qEsc[3,pInd],lty=3)
    grid()
  }
  mtext(side=2,text=yLab,outer=TRUE,line=1)
  mtext(side=1,text="Year",outer=TRUE,line=1)
}

#' Plot predicted vs observed age composition.
#'
#' @param runObj an object returned by the runJAGSmodel function.
#' @return a plot of predicted vs observed proportion age 4 escapement.
#' @export
plotAge <- function(runObj,p=0.80,newWindow=FALSE){
  qP <- (1-p)/2
  x <- getPostDraws(runObj)
  bdat <- runObj$dat$jagsDat
  odat <- runObj$dat$otherDat
  popNames <- aggregate(names(bdat$stock),list(bdat$stock),unique)$x

  # calculate expected age composition (i.e. just proportion age 4)
  pAge <- array(NA,dim=c(bdat$N,3,2))
  for(i in 1:bdat$Npops){
    pInd <- bdat$stock==i
    esc <- x$escapement[,pInd,]
    pN <- dim(esc)[2]
    tmp <- esc[,2:pN,3]/(esc[,2:pN,3] + esc[,1:(pN-1),4])
    pAge[which(pInd)[6:pN],,1] <- t(apply(tmp[,1:(pN-5)],2,quantile,prob=c(qP,0.5,1-qP)))
    pAge[which(pInd)[6:pN],,2] <- t(apply(1-tmp[,1:(pN-5)],2,quantile,prob=c(qP,0.5,1-qP)))
  }

  # calculate the credible interval for the observed data
  oAge <- array(NA,dim=c(bdat$N,3,2))
  oAge[bdat$ageInd,1,1] <- qbeta(qP,0.5+bdat$ageObs[,1],0.5+bdat$ageSampSize-bdat$ageObs[,1])
  oAge[bdat$ageInd,2,1] <- qbeta(0.5,0.5+bdat$ageObs[,1],0.5+bdat$ageSampSize-bdat$ageObs[,1])
  oAge[bdat$ageInd,3,1] <- qbeta(1-qP,0.5+bdat$ageObs[,1],0.5+bdat$ageSampSize-bdat$ageObs[,1])
  oAge[bdat$ageInd,1,2] <- qbeta(qP,0.5+bdat$ageObs[,2],0.5+bdat$ageSampSize-bdat$ageObs[,2])
  oAge[bdat$ageInd,2,2] <- qbeta(0.5,0.5+bdat$ageObs[,2],0.5+bdat$ageSampSize-bdat$ageObs[,2])
  oAge[bdat$ageInd,3,2] <- qbeta(1-qP,0.5+bdat$ageObs[,2],0.5+bdat$ageSampSize-bdat$ageObs[,2])

  # plot predicted vs observed
  if(newWindow) windows(8,8)
  par(mfrow=c(bdat$Npops,1),oma=c(3,3,1,1),mar=c(2,2,3,1))
  yrRange <- range(odat$year)
  for(i in 1:bdat$Npops){
    pInd <- bdat$stock==i
    plot(odat$year[pInd],pAge[pInd,2,1],type="l",xlim=yrRange,ylim=c(0,1),xlab="",ylab="",bty="l",main=popNames[i])
    lines(odat$year[pInd],pAge[pInd,1,1],lty=3)
    lines(odat$year[pInd],pAge[pInd,3,1],lty=3)
    segments(x0=odat$year[pInd],x1=odat$year[pInd],y0=oAge[pInd,1,1],y1=oAge[pInd,3,1],lwd=3)
    points(odat$year[pInd],oAge[pInd,2,1],pch=16)
    grid()
  }
  mtext(side=2,text="Proportion 4 year olds",outer=TRUE,line=1)
  mtext(side=1,text="Year",outer=TRUE,line=1)
}

#' Plot a parameter by year and population.
#'
#' @param runObj an object returned by the runJAGSmodel function.
#' @param param either a character string with the name of the parameter or an array with the parameter posterior.
#' @param paramName optional parameter giving the full name of the model parameter.
#' @param plotScale optional parameters describing on what scale to plot the y-axis ("identity","log","logit").
#' @return a plot of a model parameter by year and population.
#' @export
plotParam <- function(runObj,param,paramName=param,p=0.80,plotScale="identity",yRangeType="byPop",newWindow=FALSE){
  qP <- (1-p)/2
  bdat <- runObj$dat$jagsDat
  odat <- runObj$dat$otherDat
  popNames <- odat$popNames
  if(is.character(param)){
    x <- getPostDraws(runObj)
    yVals <- x[[param]]
  }else{
    yVals <- param
  }
  if(newWindow) windows(8,8)
  if(plotScale=="identity"){
    yVals <- yVals
  }else if(plotScale=="log"){
    yVals <- log(yVals)
  }else if(plotScale=="logit"){
    yVals <- logit(yVals)
  }else{
    stop("ERROR: unrecognized plotScale parameter.")
  }
  par(mfrow=c(bdat$Npops,1),oma=c(3,3,1,1),mar=c(2,2,3,1))
  qVal <- apply(yVals,2,quantile,prob=c(qP,0.5,1-qP))
  if(yRangeType=="common") yRange <- range(qVal,na.rm=TRUE)
  yrs <- odat$year
  yrRange <- range(odat$year)
  for(i in 1:bdat$Npops){
    pInd <- bdat$stock==i
    if(yRangeType=="byPop") yRange <- range(qVal[,pInd],na.rm=TRUE)
    plot(yrs[pInd],qVal[2,pInd],bty="l",xlab="", ylab="",pch=16,main=popNames[i],xlim=yrRange,ylim=yRange,type="l",lwd=2)
    lines(yrs[pInd],qVal[1,pInd],lty=3)
    lines(yrs[pInd],qVal[3,pInd],lty=3)
    grid()
  }
  mtext(side=2,text=paramName,outer=TRUE,line=1)
  mtext(side=1,text="Year",outer=TRUE,line=1)
}
