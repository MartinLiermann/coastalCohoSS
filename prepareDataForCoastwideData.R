# IN THEORY YOU ONLY NEED TO RUN THIS ONCE. BUT IT SHOULD BE RERUN WHEN CHANGE THE DATA.

# define data directory
dataDir <- "data"

# read in the OCN fish data (for all streams)
fdat2 <- read.csv(file=paste(dataDir,"fishDataBig.csv",sep="/"),stringsAsFactors=FALSE)

# create a list of the populations based on the columns
aPops <- names(fdat2)[3:23]

# change the data frame to long form
ndat <- fdat2 %>% gather(key="site",value="escapement",aPops)

# create separate wild and hatchery data frames
wdat <- ndat[ndat$hatchery=="no",c("Year","site","Harvest.Rate","escapement")]
hdat <- ndat[ndat$hatchery=="yes",c("Year","site","escapement")]
names(hdat)[names(hdat)=="escapement"] <- "escapementH"

# joint the wild and hatchery data together and create a pHOS column
tdat <- wdat %>% left_join(hdat,by=c("Year","site")) 
tdat$pHOS <- tdat$escapementH/(tdat$escapement+tdat$escapementH)
tdat$totalEscapement <- (tdat$escapement+tdat$escapementH)
tdat2 <- tdat %>% select(broodYear=Year, Site=site, escapementObs=totalEscapement, pHOS=pHOS)

# write the fish data
write.csv(tdat2,file=paste(dataDir,"fishDataOCN.csv",sep="/"),row.names=FALSE)

# create a data frame with the harvest data (not population specific)
hrdat <- wdat %>% select(broodYear=Year, site=site, HR=Harvest.Rate) %>% 
  filter(site=="Alsea") %>%
  select(broodYear=broodYear,HR=HR)

# write the harvest rate data to a file
write.csv(hrdat,file=paste(dataDir,"harvestRatesOCN.csv",sep="/"),row.names=FALSE)

# create some default habitat data and write it to a file
#  COMMENT THIS OUT WHEN WE ACTUALLY GET SOME HABITAT DATA
sites <- unique(tdat$site)
habdat <- data.frame(Site=sites,area=rep(1,length(sites)))
write.csv(habdat,file=paste(dataDir,"habitatDataOCN.csv",sep="/"),row.names=FALSE)

# create some default ocean survivals for the OCN data set 
# this needs to be data for this to work. For now just make it 1.
#  COMMENT THIS OUT WHEN IF ACTUALLY GET SOME OCEAN SURVIVAL DATA
osdat <- hrdat
osdat$OS <- rep(1,length(osdat$broodYear)) 
osdat <- osdat %>% select(broodYear,OS)
write.csv(osdat,file=paste(dataDir,"oceanSurvivalOCN.csv",sep="/"),row.names=FALSE)
