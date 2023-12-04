### Build this up backwards
## Function for sampling from AGvHD, not recent, and PR, not recent
## M paths are sampled through the multi-state model
## using sf (containing survfit elements for RFS (trans=1)
## tA is vector of length M containing times at which AGvHD was reached;
## if scalar this is interpreted as same tA for everyone
## tR is defined similarly
## c1 is the Cox model for RFS, used to adapt hazard for RFS after PR
ARsample <- function(sf,M,tA,tR,tcond,c1)
{
  # sf defined as RFS hazard at recent AGvHD and recent PR
  if (length(tA)==1) tA <- rep(tA,M)
  if (length(tR)==1) tR <- rep(tR,M)
  if (length(tcond)==1) tcond <- rep(tcond,M)
  Evs <- matrix(NA,M,2) # to contain times and events of max two steps
  tlast <- pmax(tA,tR)
  for (m in 1:M) {
    sfm <- sf
    # adapt cumulative hazard after recent periods
    sfm$haz <- diff(c(0,sfm$Haz))
    wh <- which(sfm$time > tA[m] + month)
    sfm$haz[wh] <- sfm$haz[wh] * exp(-c1$coef[["recA"]])
    wh <- which(sfm$time > tR[m] + month)
    sfm$haz[wh] <- sfm$haz[wh] * exp(-c1$coef[["recR"]])
    sfm$Haz <- cumsum(sfm$haz)
    # now sample
    crs <- mstate:::crsample(sfm, tcond=tcond[m])
    Evs[m,] <- c(crs$t,crs$trans)
  }
  return(Evs)
}
# ARsample(sf1R,M=20,tA=year*0.25-month,tR=year*0.25-month,tcond=0.25*year,c1)

Asample <- function(sf,M,tA,tcond,c1)
{
  if (length(tA)==1) tA <- rep(tA,M)
  if (length(tcond)==1) tcond <- rep(tcond,M)
  Evs <- matrix(NA,M,4) # to contain times and events of max two steps
  E1 <- matrix(NA,M,2) # for first step
  for (m in 1:M) {
    # adapt cumulative hazard after "recent" periods
    sf1 <- sf[sf$trans==1,]
    sf3 <- sf[sf$trans==3,]
    sf1$haz <- diff(c(0,sf1$Haz))
    wh <- which(sf1$time > tA[m] + month)
    sf1$haz[wh] <- sf1$haz[wh] * exp(-c1$coef[["recA"]])
    sf1$Haz <- cumsum(sf1$haz)
    sf3$haz <- diff(c(0,sf3$Haz))
    wh <- which(sf3$time > tA[m] + month)
    sf3$haz[wh] <- sf3$haz[wh] * exp(-c3$coef[["recA"]])
    sf3$Haz <- cumsum(sf3$haz)
    sfm <- rbind(sf1,sf3)    
    # now sample
    crs <- mstate:::crsample(sfm, tcond=tcond[m])
    E1[m,] <- c(crs$t,crs$trans)
  }
  Evs[,1:2] <- E1
  # Those that had platelet recovery
  whR <- which(E1[,2]==3)
  MR <- length(whR)
  if (MR>0) {
    # Hazard rate of RFS is going to be updated
    sf1R <- sf[sf$trans==1,]
    sf1R$Haz <- sf1R$Haz*exp(c1$coef[["R"]]+c1$coef[["recR"]])
    Evs[whR,3:4] <- ARsample(sf1R,M=MR,
      tA=tA[whR],tR=E1[whR,1],tcond=E1[whR,1],c1)
  }
  return(Evs)
}
# Asample(sf[sf$trans!=2,],M=20,tA=0.25*year,tcond=0.25*year,c1)

Rsample <- function(sf,M,tR,tcond,c1)
{
  if (length(tR)==1) tR <- rep(tR,M)
  if (length(tcond)==1) tcond <- rep(tcond,M)
  Evs <- matrix(NA,M,4) # to contain times and events of max two steps
  E1 <- matrix(NA,M,2) # for first step
  for (m in 1:M) {
    # adapt cumulative hazard after "recent" periods
    sf1 <- sf[sf$trans==1,]
    sf2 <- sf[sf$trans==2,]
    sf1$haz <- diff(c(0,sf1$Haz))
    wh <- which(sf1$time > tR[m] + month)
    sf1$haz[wh] <- sf1$haz[wh] * exp(-c1$coef[["recR"]])
    sf1$Haz <- cumsum(sf1$haz)
    sf2$haz <- diff(c(0,sf2$Haz))
    wh <- which(sf2$time > tR[m] + month)
    sf2$haz[wh] <- sf2$haz[wh] * exp(-c2$coef[["recR"]])
    sf2$Haz <- cumsum(sf2$haz)
    sfm <- rbind(sf1,sf2)
    # now sample
    crs <- mstate:::crsample(sf, tcond=tcond[m])
    E1[m,] <- c(crs$t,crs$trans)
  }
  Evs[,1:2] <- E1
  # Those that had AGvHD
  whA <- which(E1[,2]==2)
  MA <- length(whA)
  if (MA>0) {
    # Hazard rate of RFS is going to be updated
    sf1A <- sf[sf$trans==1,]
    sf1A$Haz <- sf1A$Haz*exp(c1$coef[["A"]]+c1$coef[["recA"]])
    Evs[whA,3:4] <- ARsample(sf1A,M=MA,
      tA=E1[whA,1],tR=tR[whA],tcond=E1[whA,1],c1)
  }
  return(Evs)
}
# Rsample(sf[sf$trans!=3,],M=20,tR=year*0.15-month,tcond=0.15,c1)

Msample <- function(sf,M,tcond,c1,c2,c3)
{
  if (length(tcond)==1) tcond <- rep(tcond,M)
  Evs <- matrix(NA,M,6) # to contain times and events of max two steps
  E1 <- matrix(NA,M,2) # for first step
  for (m in 1:M) {
    crs <- mstate:::crsample(sf,tcond=tcond[m])
    E1[m,] <- c(crs$t,crs$trans)
  }
  Evs[,1:2] <- E1
  ### First event generated
  # is.na(E1[,2]) have no event at all, no further action required
  # E1[,2]==1 have relapse or death as event, no further action required
  # E1[,2]==2 have AGvHD as first intermediate event
  # these will sampled further with recAsample
  whA <- which(E1[,2]==2)
  MA <- length(whA)
  if (MA>0) {
    sf1 <- sf[sf$trans==1,]
    sf1$Haz <- sf1$Haz * exp(c1$coef[["A"]]+c1$coef[["recA"]])
    sf3 <- sf[sf$trans==3,]
    sf3$Haz <- sf3$Haz * exp(c3$coef[["A"]]+c3$coef[["recA"]])
    sfA <- rbind(sf1,sf3)
    Evs[whA,3:6] <- Asample(sfA,MA,E1[whA,1],E1[whA,1],c1)
  }
  # E1[,2]==3 have PR as first intermediate event
  # these will sampled further with recRsample
  whR <- which(E1[,2]==3)
  MR <- length(whR)
  if (MR>0) {
    sf1 <- sf[sf$trans==1,]
    sf1$Haz <- sf1$Haz * exp(c1$coef[["R"]]+c1$coef[["recR"]])
    sf2 <- sf[sf$trans==2,]
    sf2$Haz <- sf2$Haz * exp(c2$coef[["R"]]+c2$coef[["recR"]])
    sfR <- rbind(sf1,sf2)
    Evs[whR,3:6] <- Rsample(sfR,MR,E1[whR,1],E1[whR,1],c1)
  }
  return(Evs)
}

extractRFS <- function(Events,t)
{
  M <- nrow(Events)
  m <- ncol(Events)/2
  cols <- (1:m)*2 - 1
  whc <- apply(!is.na(Events[,cols,drop=FALSE]),1,sum)
  RFS <- rep(NA,M)
  RFS[whc==1] <- Events[whc==1,1]
  RFS[whc==2] <- Events[whc==2,3]
  RFS[whc==3] <- Events[whc==3,5]
  RFS <- data.frame(time=sort(RFS,decreasing=TRUE),prob=rev((1:M)/M))
  RFS <- RFS[!duplicated(RFS$time),]
  RFS <- data.frame(time=rev(RFS$time),prob=rev(RFS$prob))
  RFS <- RFS[!is.infinite(RFS$time),]
  return(evalstep(RFS$time,RFS$prob,t))
}
