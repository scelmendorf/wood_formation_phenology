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
lat<- 47
fstar <- 200
cstar <- 110
fstaradjforchill <- 3 # how much more GDD to you need based on diff from cstar at end of daystostart
basetemp <- 6 #mean annual temp
sigma <- 4

# 2. Put together the seasonal temps, varying fstar (increases when chill is low) and calculate sensitivities

set.seed(113)

df <- data.frame(lat = numeric(), site = numeric(), year = numeric(), 
                 mat = numeric(), gdd = numeric(), photo = numeric(),
                 phen = numeric())  

  for (j in 1:sitez){
    yearly_expected_temp <- rep(basetemp, yearz)
    daily_temp <- sapply(yearly_expected_temp, function(x) c(rnorm(daysperseason, 0, sigma),
                                                             rnorm(daysperinterseason, 2, sigma), rnorm(daysperinterseason, 4, sigma),
                                                             rnorm(daysperseason, 6, sigma)))#daily mean temperature from winter to spring (warms from mean = 0 to mean = 6)
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
    
    #plot relationships       
    figname<-paste("figures/onelat",lat,yearz, "years","site",j,".pdf",sep = "_")
    pdf(figname,height = 6,width = 8)
    par(mfrow=c(1,3))
    plot(yearly_tempall,phen_doy, pch=16, bty = "l",xlab = "Temperature", ylab = "Day of Year", main = paste("site",j, sep=""))
    m<-lm(phen_doy~yearly_tempall)
    if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2)}
    
    plot(gddsum_trunc, phen_doy,pch=16,col="darkgreen", bty="l",xlab="GDD",ylab="Day of Year", main = paste("site",j, sep=""))
    m<-lm(phen_doy~gddsum_trunc)
    if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkgreen")}
    
    plot(dl, phen_doy,pch=16,col="darkblue", bty="l",xlab="Daylength (hrs)",ylab="Day of Year", main = paste("site",j, sep=""))
    m<-lm(phen_doy~dl)
    if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkblue")}
    
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
latz<-seq(from=20, to = 70, by = 10)
fstar <- 200
cstar <- 110
fstaradjforchill <- 3 # how much more GDD to you need based on diff from cstar at end of daystostart
basetemp <- 6 #mean annual temp depends on latitude. ???0.87 °C / per degree latitude (Wang et al journal of mtn science)
sigma <- 4


# 2. Put together the seasonal temps, varying fstar (increases when chill is low) and calculate sensitivities

 set.seed(113)

 df <- data.frame(lat = numeric(), site = numeric(), year = numeric(), 
                    mat = numeric(), gdd = numeric(), photo = numeric(),
                    phen = numeric())  


#figname<-paste("figures/plot_site",site2plot,".pdf",sep="_")
#pdf(figname, width = 12, height = 3))
for (l in 1:length(latz)){
  for (j in 1:sitez){
     yearly_expected_temp <- rep(basetemp, yearz)
     daily_temp <- sapply(yearly_expected_temp, function(x) c(rnorm(daysperseason, 0, sigma),
           rnorm(daysperinterseason, 2, sigma), rnorm(daysperinterseason, 4, sigma),
           rnorm(daysperseason, 6, sigma)))#daily mean temperature from winter to spring (warms from mean = 0 to mean = 6)
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
      windows()
      par(mfrow=c(1,3))
      plot(yearly_tempall,phen_doy, pch=16, bty = "l",xlab = "Temperature", ylab = "Day of Year", main = paste("site",j, sep=""))
      m<-lm(phen_doy~yearly_tempall)
      if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2)}
      
      plot(gddsum_trunc, phen_doy,pch=16,col="darkgreen", bty="l",xlab="GDD",ylab="Day of Year", main = paste("site",j, sep=""))
      m<-lm(phen_doy~gddsum_trunc)
      if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkgreen")}
      
      plot(dl, phen_doy,pch=16,col="darkblue", bty="l",xlab="Daylength (hrs)",ylab="Day of Year", main = paste("site",j, sep=""))
      m<-lm(phen_doy~dl)
      if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkblue")}
      
      
  
  #save all the stuff we care about   
  dfadd <- data.frame(lat = latz[l], site=j, year = seq(1:yearz), 
                      mat =  yearly_tempall, gdd=gddsum_trunc, photo = dl,
                      phen=phen_doy)    
  df <- rbind(df, dfadd)
  }#j
  
}#l

 windows()
 par(mfrow=c(2,3))
 plot(df$mat,df$phen, pch=16, bty = "l",xlab = "Temperature", ylab = "Day of Year")
 m<-lm(df$phen~df$mat)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2)}
 
 plot(df$gdd, df$phen,pch=16,col="darkgreen", bty="l",xlab="GDD",ylab="Day of Year")
 m<-lm(df$phen~df$gdd)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkgreen")}
 
 plot(df$photo, df$phen,pch=16,col="darkblue", bty="l",xlab="Daylength (hrs)",ylab="Day of Year")
 m<-lm(phen_doy~dl)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkblue")}
 
 plot(df$photo, df$mat,pch=16,col="gray", bty="l",xlab="Daylength (hrs)",ylab="Mean Annual Temperature")
 m<-lm(phen_doy~dl)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="gray")}
   
   
 