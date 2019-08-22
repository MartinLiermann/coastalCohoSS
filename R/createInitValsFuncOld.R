#' create function that generates initial values for the Grand Ronde chinook model
#'
#' @param bdat the model data.
#' @return a function which generates a list of initial values.
#' @examples
#' initFunc <- createInitValsFunc(bdat)
#'
#' @export
#
####################################################
# WHAT I NEED TO DO NEXT
####################################################
createInitValsFuncOld <- function(dat, includePlots=FALSE){
  bdat <- dat$jagsDat
  popNames <- tapply(names(bdat$stock),list(bdat$stock),unique)
  function(){
    logit <- function(x) log(x/(1-x))
    invLogit <- function(x) 1/(1+exp(-x))

    ##### expand the observation data to include all combinations
    #####   of years and populations so it matches the parameters
    edat <- expandAllData(bdat)

    #####  spawners to trap model (prepare values) #####
    spawners <- edat$escapementFemObs
    smolt <- edat$smoltObs
    habVar <- bdat$habVar

    ##### estimate the SR params #####
    srF <- function(p,ss) ss/(1/p[1] + ss/p[2])
    fitBH <- function(ss,rr){
      ssq <- function(pp,ss,rr){
        p <- exp(pp)
        sum((rr-srF(p,ss))^2)
      }
      bFit <- nlm(ssq,p=log(c(median(rr/ss),quantile(rr,prob=0.8))),ss,rr)
      bFit
    }
    habDat <- habVar[edat$stock]
    escDat <- spawners/(1-edat$pHOS)
    escDat <- pmax(1,escDat)
    srDat <- data.frame(S=escDat/habDat, R=smolt/habDat, stock=edat$stock)
    srDat <- srDat[!is.na(srDat$S) & !is.na(srDat$R),] # remove NAs
    srParams <- array(NA,dim=c(bdat$Npops,2))
    if(includePlots){
      windows(10,10)
      par(mfrow=c(3,2))
      xx <- 0:round(max(srDat$S))
    }
    xLim <- c(0,max(srDat$S))
    yLim <- c(0,max(srDat$R))
    for(st in 1:bdat$Npops){
      sInd <- srDat$stock==st
      r1 <- fitBH(srDat$S[sInd],srDat$R[sInd])
      srParams[st,] <- exp(r1$estimate)
      if(includePlots){
        plot(srDat$S[sInd],srDat$R[sInd],pch=16,col="black",xlab="S/PEU",ylab="Parr/PEU",bty="l",main=popNames[st],xlim=xLim,ylim=yLim)
        lines(xx,srF(exp(r1$estimate),ss=xx),lwd=2,col=rgb(0.2,0.8,0.2,0.5))
        grid()
      }
    }

    # now calculate the model specific parameters
    compTau <- function(x) 1/var(x,na.rm=T)
    prod <- srParams[,1]
    cap <- srParams[,2]
    capSlope <- median(cap/bdat$habVar)
    escObs <- edat$escapementFemObs
    esc <- rep(NA,length(escObs))
    for(i in 1:bdat$Npops){
      pInd <- which(edat$stock==i & edat$year %in% (1:3))
      esc[pInd] <- escObs[pInd]
    }

    pInd <- bdat$stock
    initList <- list(
      prod = prod,
      cap = cap,
      capSlope = capSlope,
      yearEffectTmp = rep(0,bdat$Nyears)
      #spawnersFem = esc
   )
    initList
  }
}
