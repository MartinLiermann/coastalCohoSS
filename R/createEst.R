# Function to create rough estimates of parameters based on data
createEst <- function(dataDir){
  rdat <- read.csv(file=paste(dataDir,"together.csv",sep=""),stringsAsFactors=FALSE)
  hdat <- read.csv(file=paste(dataDir,"habitatData.csv",sep=""),stringsAsFactors=FALSE)
  habDat1 <- hdat$poolPer10k
  names(habDat1) <- hdat$population
  habDat <- habDat1[rdat$population]
  pops <-hdat$population
  parr <- (rdat$spring_survToLGR*rdat$spr_mig + rdat$fall_survToLGR*rdat$fall_mig)/
    rdat$summer_survToLGR

  pOut <- rdat$fall_mig/parr
  winterSurv <- pmin(0.999,rdat$spr_mig / (parr*(1-pOut)))
  moveSurvRes <- rdat$spring_survToLGR
  moveSurvMig <- pmin(0.999,rdat$fall_survToLGR*rdat$fall_mig / (parr * pOut * winterSurv))

  # we have two different ways of estimating winter survival (what we did above, and dividing the observed winter parr to dam survival by moveSurvRes)
  winterSurv2 <- rdat$winter_survToLGR/rdat$spring_survToLGR
  survMax <- max(c(winterSurv,winterSurv2),na.rm=T)

  # from Tom's stuff
  Spr_LGR <- rdat$spring_survToLGR*rdat$spr_mig
  win_parr <- Spr_LGR/rdat$winter_survToLGR

  pseudoparr <- win_parr + rdat$fall_mig

  srDat <- data.frame(
    S=rdat$spawners45*rdat$Agewt/habDat,
    parr=parr/habDat,
    ws=winterSurv,
    ws2=winterSurv2,
    ws3=win_parr/parr,
    win_parr=win_parr,
    pseudoparr=pseudoparr,
    pop=rdat$population)

  plotVars <- function(srDat,x,y,xLab=NULL,yLab=NULL,xLim=NULL,yLim=NUll){
    if(is.null(xLim)) xLim <- range(srDat[,x],na.rm=TRUE)
    if(is.null(yLim)) yLim <- range(srDat[,y],na.rm=TRUE)
    if(is.null(xLab)) xLab <- x
    if(is.null(yLab)) yLab <- y
    windows(10,10)
    par(mfrow=c(2,2))
    pops <- unique(srDat$pop)
    for(pop in pops){
      pDat <- srDat[srDat$pop==pop,]
      plot(pDat[,x],pDat[,y],xlab=xLab,ylab=yLab,bty="l",main=pop,
           pch=16,col=rgb(0,0,0,0.5),xlim=xLim,ylim=yLim)
    }
  }

  plotVars(srDat,"S","parr",xLab="Spawner/PEU",yLab="Parr/PEU",xLim=c(0,100),yLim=c(0,20000))
  plotVars(srDat,"parr","ws",xLab="Parr/PEU",yLab="Winter Survival",xLim=NULL,yLim=c(0,1))
#  plotVars(srDat,"parr","ws2",xLab="Parr/PEU",yLab="Winter Survival 2",xLim=NULL,yLim=c(0,1))
#  plotVars(srDat,"parr","ws3",xLab="Parr/PEU",yLab="Winter Survival 3",xLim=NULL,yLim=c(0,1))

}
