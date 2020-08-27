## Started 24 Aug 2929 ##
## By Ailene ##
## Modified from Lizzie's declining sensitivity code
## This code requires higher forcing (aka GDD, aka fstar) for leafout when chilling is low #
## It does not vary when GDD accumulates though, it always just starts accumulating on a certain day (daystostart) #
## Then looks at how well photoperiod does at predicting day of year of the phenological event
## Note: photoperiod is not used to simulate phenology doy data at all 
## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

#load librarys
library(insol)
library(scales)

# Setting working directory. 
if(length(grep("ailene", getwd()))>0) { 
  setwd("C:/Users/ailene.ettinger/Documents/GitHub/wood_formation_phenology")
} 


set.seed(113)

## Simulate daily temperature and phenology data
# 1. First, for one latitude only (with multiple sites and years): set up years, seasons (winter and spring), days per season, temperatures, required GDD (fstar), required chill (cstar) and how much higher fstar is when cstar is not met
daysperseason <- 100
daysperinterseason <- 25
daystostart <- daysperseason+daysperinterseason # this defines the break between 'winter' and 'spring,' functionally we accumulate chill only in 1:daystostart and GDD only in daystostart:end(df)
yearz <- 20
sitez <- 10
lat<- 65
fstar <- 200
cstar <- 110
fstaradjforchill <- 3 # how much more GDD to you need based on diff from cstar at end of daystostart
basetemp <- 6 #mean annual temp
sigma <- 8

# 2. Put together the seasonal temps, varying fstar (increases when chill is low) and calculate sensitivities

set.seed(113)

df <- data.frame(lat = numeric(), site = numeric(), year = numeric(), 
                 mat = numeric(), gdd = numeric(), photo = numeric(),
                 phen = numeric())  

  for (j in 1:sitez){
    yearly_expected_temp <- rep(basetemp, yearz)
    daily_temp <- sapply(yearly_expected_temp, function(x) c(rnorm(daysperseason, basetemp - 6, sigma),
                                                             rnorm(daysperinterseason, basetemp - 4, sigma), rnorm(daysperinterseason, basetemp - 2, sigma),
                                                             rnorm(daysperseason, basetemp, sigma)))#daily mean temperature from winter to spring (warms from mean = 0 to mean = 6)
    chill <- daily_temp #set up chilling
    chill[(chill)<0] <- 0
    chill[(chill)>5] <- 0
    gdd <- daily_temp #set up gdd
    gdd[(gdd)<0] <- 0
    gddreq <- c()
    phen_dat <- c()
    for (k in 1:ncol(chill)){
      chillsum <- sapply(1:ncol(chill), function(x) (cumsum(chill[,x])))
      gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[daystostart:nrow(gdd),x])))
      if (chillsum[daystostart,k]>cstar) {
        gddreq[k] <- fstar
      } else {
        gddreq[k] <- fstar + (cstar-chillsum[daystostart,k])*fstaradjforchill
      }
      phen_dat[k] <- min(which(gddsum[,k] > gddreq[k]))
      meanchill <- mean(chillsum[daystostart,])#why taking mean here? mean across 30 years?
      meanfstar <- mean(gddreq)
      chillmet<-length(which(chillsum[daystostart,]>cstar))/yearz
      meangddsum<- mean(gddsum)
    }#k
    phen_doy<-phen_dat+daystostart#convert to doy of year
    yearly_tempall <- colMeans(daily_temp)
    yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:phen_dat[x], x]))
    gddsum_trunc <- sapply(1:ncol(gddsum), function(x) gddsum[phen_dat[x],x])
    dl<-as.data.frame(daylength(lat,-120,phen_doy,8))$daylen#longitude, timezome for Seattle
    all_doy<-seq(from =daystostart, to =(nrow(gddsum)+daystostart), by = 1)
    all_dl<-as.data.frame(daylength(lat,-120,all_doy,8))$daylen#
    #plot relationships       
    figname<-paste("figures/onelat",lat,yearz, "years","site",j,".pdf",sep = "_")
    pdf(figname,height = 4,width = 12)
    par(mfrow=c(1,3))
    plot(yearly_tempall,phen_doy, pch=16, bty = "l",xlab = "Temperature", ylab = "Day of Year",ylim = range(all_doy), main = paste("site",j, sep=""))
    m<-lm(phen_doy~yearly_tempall)
    if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2)}
    
    plot(gddsum_trunc, phen_doy,pch=16,col="darkgreen", bty="l",xlab="GDD",ylab="Day of Year", ylim = range(all_doy),main = paste("site",j, sep=""))
    m<-lm(phen_doy~gddsum_trunc)
    if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkgreen")}
    #lines(gddsum, all_doy)
    plot(dl, phen_doy,pch=16,col="darkblue", bty="l",xlab="Daylength (hrs)",ylab="Day of Year", ylim = range(all_doy),xlim = range(all_dl),main = paste("site",j, sep=""))
    m<-lm(phen_doy~dl)
    if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkblue")}
    lines(all_dl, all_doy, col = alpha("darkblue", alpha = .2))
    
    dev.off()
    
    #save all the stuff we care about   
    dfadd <- data.frame(lat = lat, site=j, year = seq(1:yearz), 
                        mat =  yearly_tempall, gdd=gddsum_trunc, photo = dl,
                        phen=phen_doy)    
    df <- rbind(df, dfadd)
  }#j


#2. Now for multiple sites, at different latitudes  

daysperseason <- 100
daysperinterseason <- 25
daystostart <- daysperseason+daysperinterseason # this defines the break between 'winter' and 'spring,' functionally we accumulate chill only in 1:daystostart and GDD only in daystostart:end(df)
yearz <- 20
sitez <- 10
latz<-seq(from=40, to = 70, by = 5)
fstar <- 180
# if we  allow fstar to vary with latitude 
#latfstareff<-#not sure what to put here- will look at some refs> GDD_threshold_(persite) ~ intercept + B*site_latitude
cstar <- 110
fstaradjforchill <- 3 # how much more GDD to you need based on diff from cstar at end of daystostart
basetemp <- 25#mean annual temp (this is warmer than above because its an intercept now, for a lm with lat as a predictor for temp
lattempeff<- ???0.87# °C / per degree latitude (Wang et al journal of mtn science)
sigma <- 8


# 2. Put together the seasonal temps, varying fstar (increases when chill is low) and calculate sensitivities

 set.seed(113)

 dflat <- data.frame(lat = numeric(), site = numeric(), year = numeric(), 
                    mat = numeric(), gdd = numeric(), photo = numeric(),
                    phen = numeric())  


#figname<-paste("figures/plot_site",site2plot,".pdf",sep="_")
#pdf(figname, width = 12, height = 3))
for (l in 1:length(latz)){
  for (j in 1:sitez){
    meantemp<- basetemp+(latz[l]-latz[1])*lattempeff
    yearly_expected_temp <- rep(meantemp, yearz)
     daily_temp <- sapply(yearly_expected_temp, function(x) c(rnorm(daysperseason, meantemp-6, sigma),
           rnorm(daysperinterseason, meantemp - 4, sigma), rnorm(daysperinterseason, meantemp - 2, sigma),
           rnorm(daysperseason, meantemp, sigma)))#daily mean temperature from winter to spring (warms from mean = 0 (or basetemp - 6) to mean = 6 (or whatever basetemp is))
     chill <- daily_temp #set up chilling
       chill[(chill)<0] <- 0
       chill[(chill)>5] <- 0
      gdd <- daily_temp #set up gdd
       gdd[(gdd)<0] <- 0
       gddreq <- c()
       phen_dat <- c()
       for (k in 1:ncol(chill)){
           chillsum <- sapply(1:ncol(chill), function(x) (cumsum(chill[,x])))
           gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[daystostart:nrow(gdd),x])))
           if (chillsum[daystostart,k]>cstar) {
           gddreq[k] <- fstar
           } else {
           gddreq[k] <- fstar + (cstar-chillsum[daystostart,k])*fstaradjforchill
           }
          phen_dat[k] <- min(which(gddsum[,k] > gddreq[k]))
           meanchill <- mean(chillsum[daystostart,])#why taking mean here? mean across 30 years?
           meanfstar <- mean(gddreq)
           chillmet<-length(which(chillsum[daystostart,]>cstar))/yearz
           meangddsum<- mean(gddsum)
           }#k
    phen_doy<-phen_dat+daystostart#convert to doy of year
    yearly_tempall <- colMeans(daily_temp)
    yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:phen_dat[x], x]))
    gddsum_trunc <- sapply(1:ncol(gddsum), function(x) gddsum[phen_dat[x],x])
    dl<-as.data.frame(daylength(latz[l],-120,phen_doy,8))$daylen#longitude, timezome for Seattle
    
    #plot relationships       
      #windows()
      # par(mfrow=c(1,3))
      # plot(yearly_tempall,phen_doy, pch=16, bty = "l",xlab = "Temperature", ylab = "Day of Year", main = paste("site",j, sep=""))
      # m<-lm(phen_doy~yearly_tempall)
      # if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2)}
      # 
      # plot(gddsum_trunc, phen_doy,pch=16,col="darkgreen", bty="l",xlab="GDD",ylab="Day of Year", main = paste("site",j, sep=""))
      # m<-lm(phen_doy~gddsum_trunc)
      # if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkgreen")}
      # 
      # plot(dl, phen_doy,pch=16,col="darkblue", bty="l",xlab="Daylength (hrs)",ylab="Day of Year", main = paste("site",j, sep=""))
      # m<-lm(phen_doy~dl)
      # if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkblue")}
      # 
      
  #save all the stuff we care about   
  dflatadd <- data.frame(lat = latz[l], site=j, year = seq(1:yearz), 
                      mat =  yearly_tempall, gdd=gddsum_trunc, photo = dl,
                      phen=phen_doy)    
  dflat <- rbind(dflat, dflatadd)
  }#j
  
}#l

 #Try making some plots like the authors
 figname<-"figures/alllatsallyears.pdf"
 pdf(figname,height = 6,width = 10)
 
 par(mfrow=c(2,3))
 plot(dflat$mat,dflat$phen, pch=16, bty = "l",xlab = "Temperature", ylab = "Day of Year")
 m<-lm(dflat$phen~dflat$mat)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2)}
 
 plot(dflat$gdd, dflat$phen,pch=16,col="darkgreen", bty="l",xlab="GDD",ylab="Day of Year")
 m<-lm(dflat$phen~dflat$gdd)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkgreen")}
 
 plot(dflat$photo, dflat$phen,pch=16,col="darkblue", bty="l",xlab="Daylength (hrs)",ylab="Day of Year")
 m<-lm(dflat$phen~dflat$photo)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkblue")}
 
 plot(dflat$photo, dflat$mat,pch=16,col="gray", bty="l",xlab="Daylength (hrs)",ylab="Mean Annual Temperature")
 m<-lm(dflat$mat~dflat$photo)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="gray")}
   
 plot(dflat$lat, dflat$mat,pch=16,col="gray", bty="l",xlab="Latitude (degrees)",ylab="Mean Annual Temperature")
 m<-lm(dflat$mat~dflat$lat)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="gray")}
 plot(dflat$lat, dflat$photo,pch=16,col="gray", bty="l",xlab="Latitude (degrees)",ylab="Daylength (hrs)")
 m<-lm(dflat$mat~dflat$photo)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="gray")}
dev.off()
#check that temp is changing with latitude:
aggregate(dflat$mat, by=list(dflat$lat), mean)