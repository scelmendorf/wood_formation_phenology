## Started 24 Aug 2929 ##
## By Ailene ##
## Modified from Lizzie's declining sensitivity code
## This code requires a threshold amount of forcing (aka GDD, aka fstar) #
## It does not vary when GDD accumulates or what the threshold is (e.g. with latitude), 
# it always just starts accumulating on a certain day (Jan 1, as in Huang et al 2020) #
## Then looks at how well photoperiod does at predicting day of year of the phenological event
## Note: photoperiod is not used to simulate phenology doy data at all 
## housekeeping

rm(list=ls()) 
options(stringsAsFactors = FALSE)

#load librarys
library(insol)
library(scales)
library(RColorBrewer)

# Setting working directory. 
if(length(grep("ailene", getwd()))>0) { 
  setwd("C:/Users/ailene.ettinger/Documents/GitHub/wood_formation_phenology")
} 


## Simulate daily temperature and phenology data

daysperseason <- 60
daysperinterseason <- 30
yearz <- 20
sitez <- 10
fstar <- 400 

lat<- 45
basetemp <- 15
lattempeff<- ???0.87
meantemp<- basetemp+(lat-lat)*lattempeff
sigma <- 8

# 2. Simulate winter-spring daily temperature, calculate forcing, and phenology for one site

df <- data.frame(lat = numeric(), site = numeric(), year = numeric(), 
                 mat = numeric(), fu = numeric(), photo = numeric(),
                 phen = numeric(), phenrand = numeric())  

for (j in 1:sitez){
    yearly_expected_temp <- rep(basetemp, yearz)
    daily_temp <- sapply(yearly_expected_temp, function(x) c(rnorm(daysperseason, meantemp-6, sigma),
                                                             rnorm(daysperinterseason, meantemp - 4, sigma), 
                                                             rnorm(daysperinterseason, meantemp - 2, sigma),
                                                             rnorm(daysperseason, meantemp, sigma)))#daily mean temperature from winter to spring (warms from mean = 0 (or basetemp - 6) to mean = 6 (or whatever basetemp is))
    
    #plot(seq(1:dim(daily_temp)[1]),daily_temp[,1], type = "l")
         
    force = ifelse(daily_temp>5, 28.4/(1+exp(-.185*(daily_temp-18.4))),0)
   
     fureq <- c()
    phen_dat <- c()
    for (k in 1:yearz){
      forcesum <- sapply(1:ncol(force), function(x) (cumsum(force[1:nrow(force),x])))
      fureq[k] <- fstar
      phen_dat[k] <- min(which(forcesum[,k] > fureq[k]))
      meanfstar <- mean(fureq)
      meanforcesum<- mean(forcesum)
    }#k
    phen_doy<-phen_dat#convert to doy of year
  
    yearly_tempall <- colMeans(daily_temp)
    yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:phen_dat[x], x]))
    forcesum_trunc <- sapply(1:ncol(forcesum), function(x) forcesum[phen_dat[x],x])
    dl<-as.data.frame(daylength(lat,-120,phen_doy,8))$daylen#longitude, timezome for Seattle
    all_doy<-seq(from =1, to =(nrow(forcesum)), by = 1)
    all_dl<-as.data.frame(daylength(lat,-120,all_doy,8))$daylen#
    
    #instead of using phenology that depends on  GDD, what if we just make it random?
    phen_doy_rand<-rnorm(yearz,mean(phen_doy),100)#
    dl_rand<-as.data.frame(daylength(lat,-120,phen_doy_rand,8))$daylen#longitude, timezome for Seattle

    if(j == 1){figname<-paste("plots/onelat",lat,yearz, "years","site",j,".png",sep = "_")
      png(figname,height = 400,width = 1400,res = 100)
      par(mfrow=c(1,3))
      
      plot(forcesum_trunc, phen_doy,pch=16,col="darkgreen", 
           xlim = range(forcesum), ylim = range(all_doy),bty="l",
           xlab="Forcing",ylab="Onset date of wood formation (DOY)",
           cex.lab = 2, cex.axis = 2)
      for (l in 1:yearz){
      lines(forcesum[,l],all_doy, col = "darkgray", lwd = 1.5)
      }
      abline(v= 400, lty = 2, lwd = 2, col = "darkgreen")
      points(forcesum_trunc, phen_doy,pch=16,cex = 1.5, col = "darkgreen")
      
      text(fstar+100,max(all_doy),"FU", cex = 2)
      
      mtext("A)", line = 2, adj= 0)
      
      plot(yearly_tempall, phen_doy,pch=16,col="darkgray", bty="l",
           xlab="Mean Annual Temperature (C))",ylab="Onset date of wood formation (DOY)",
           cex.lab = 2, cex.axis = 2)
      m<-lm(phen_doy~yearly_tempall)
      if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="blue")}
      mtext(paste("R2 = ",round(summary(m)$r.squared, digits = 2), sep = ""), side = 1,adj = 1, line = -2)
      mtext("B)", line = 2, adj= 0)
      
      plot(dl, phen_doy,pch=16,col="darkgray", bty="l",
           xlab="Photoperiod (hrs)",ylab="Onset date of wood formation (DOY)",
           cex.lab = 2, cex.axis = 2)
      m<-lm(phen_doy~dl)
      if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="blue")}
      mtext(paste("R2 = ",round(summary(m)$r.squared, digits = 2), sep = ""), side = 1,adj = 1, line = -2)
      mtext("C)", line = 2, adj= 0)
      
      # plot(dl_rand, phen_doy_rand,pch=16,col="darkgray", bty="l",
      #      xlab="Daylength (hrs)",ylab="Onset date of wood formation (DOY)",ylim = range(all_doy),xlim = range(all_dl),
      #      cex.lab = 2, cex.axis = 2)
      # lines(all_dl, all_doy, col = alpha("darkgray", alpha = .4), lwd = 2)
      # m<-lm(phen_doy_rand~dl_rand)
      # clip(min(dl_rand), max(dl_rand),min(phen_doy_rand), max(phen_doy_rand))
      # if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="blue")}
      # mtext("D)", line = 2, adj= 0)

    dev.off()
    }
    
    #save all the stuff we care about   
    dfadd <- data.frame(lat = lat, site=j, year = seq(1:yearz), 
                        mat =  yearly_tempall, fu=forcesum_trunc, photo = dl,
                        phen=phen_doy, phenrand = phen_doy_rand)    
    df <- rbind(df, dfadd)
  }#j


#2. Now for multiple sites, at different latitudes  

daysperseason <- 60
daysperinterseason <- 30
yearz <- 20
sitez <- 10
basetemp <- 18
lattempeff<- -0.87

latz<-seq(from=40, to = 60, by = 5)
fstar <- 172
lattempeff<- -0.5# °C / per degree latitude (Wang et al journal of mtn science)
sigma <- 8

figname<-paste("plots/onelatalllats.png")
png(figname,height = 800,width = 1400, res = 100)
par(mfrow=c(2,3), mar=c(5,5,4,2))

 dflat <- data.frame(lat = numeric(), site = numeric(), year = numeric(), 
                    mat = numeric(), fu = numeric(), photo = numeric(),
                    phen = numeric())  
  
for (l in 1:length(latz)){
  for (j in 1:sitez){
    meantemp<- basetemp+(latz[l]-latz[1])*lattempeff
    yearly_expected_temp <- rep(meantemp, yearz)
     daily_temp <- sapply(yearly_expected_temp, function(x) c(rnorm(daysperseason, meantemp-6, sigma),
           rnorm(daysperinterseason, meantemp - 4, sigma), rnorm(daysperinterseason, meantemp - 2, sigma),
           rnorm(daysperseason+20, meantemp, sigma)))#daily mean temperature from winter to spring (warms from mean = 0 (or basetemp - 6) to mean = 6 (or whatever basetemp is))
      
      
       #plot(seq(1:dim(daily_temp)[1]),daily_temp[,1], type = "l")
       
       force = ifelse(daily_temp>5, 28.4/(1+exp(-.185*(daily_temp-18.4))),0)
       
       fureq <- c()
       phen_dat <- c()
       for (k in 1:yearz){
         forcesum <- sapply(1:ncol(force), function(x) (cumsum(force[1:nrow(force),x])))
         fureq[k] <- fstar
         #if(length(which(forcesum[,k] > fureq[k]))==0){phen_dat[k]<-dim(daily_temp)[1]}
         #if (length(which(forcesum[,k] > fureq[k]))>0){
           phen_dat[k] <- min(which(forcesum[,k] > fureq[k]))#}
         meanfstar <- mean(fureq)
         meanforcesum<- mean(forcesum)
       }#k
       phen_doy<-phen_dat+30#convert to doy of year

    yearly_tempall <- colMeans(daily_temp)
    
    yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:phen_dat[x], x]))
    forcesum_trunc <- sapply(1:ncol(forcesum), function(x) forcesum[phen_dat[x],x])
    all_doy<-seq(from = 30, to = dim(daily_temp)[1]+29, by = 1)
    dl<-as.data.frame(daylength(latz[l],-120,phen_doy,8))$daylen#longitude, timezome for Seattle
      if( latz[l]==45 & j == 1){
        plot(forcesum_trunc, phen_doy,pch=16,col="darkgreen", 
         xlim = range(forcesum), ylim = range(all_doy),bty="l",
         xlab="Forcing",ylab="Onset date of wood formation (DOY)",
         cex.lab = 2, cex.axis = 2)
        for (y in 1:yearz){
      lines(forcesum[,y],all_doy, col = "darkgray", lwd = 1.5)}
      abline(v= fstar, lty = 2, lwd = 2, col = "darkgreen")
      points(forcesum_trunc, phen_doy,pch=16,cex = 1.5, col = "darkgreen")
    
      text(fstar,max(all_doy),"FU", cex = 2)
    
      mtext("A)", line = 2, adj= 0)
    
    plot(yearly_tempall, phen_doy,pch=16,col="darkgray", bty="l",
         xlab="Mean Annual Temperature (C))",ylab="Onset date of wood formation (DOY)",
         cex.lab = 2, cex.axis = 2)
        m<-lm(phen_doy~yearly_tempall)
      if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="blue")}
      mtext(paste("R2 = ",round(summary(m)$r.squared, digits = 2), sep = ""), side = 1,adj = 1, line = -2)
      mtext("B)", line = 2, adj= 0)
    
    plot(dl, phen_doy,pch=16,col="darkgray", bty="l",
         xlab="Photoperiod (hrs)",ylab="Onset date of wood formation (DOY)",
         cex.lab = 2, cex.axis = 2)
    m<-lm(phen_doy~dl)
    if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="blue")}
    mtext(paste("R2 = ",round(summary(m)$r.squared, digits = 2), sep = ""), side = 1,adj = 1, line = -2)
    mtext("C)", line = 2, adj= 0)
    
     # plot(dl_rand, phen_doy_rand,pch=16,col="darkgray", bty="l",
     #      xlab="Photoperiod (hrs)",ylab="Onset date of wood formation (DOY)",ylim = range(all_doy),xlim = range(all_dl),
     #      cex.lab = 2, cex.axis = 2)
     # lines(all_dl, all_doy, col = alpha("darkgray", alpha = .4), lwd = 2)
     # m<-lm(phen_doy_rand~dl_rand)
     # clip(min(dl_rand), max(dl_rand),min(phen_doy_rand), max(phen_doy_rand))
     # if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="blue")}
     # mtext("D)", line = 2, adj= 0)

    
        }
  
      
  #save all the stuff we care about   
  dflatadd <- data.frame(lat = latz[l], site=j, year = seq(1:yearz), 
                      mat =  yearly_tempall, fu=forcesum_trunc, photo = dl,
                      phen=phen_doy)    
  dflat <- rbind(dflat, dflatadd)
  }#j
  
}#l

 plot(dflat$fu, dflat$phen,pch=16,col="darkgreen", bty="l",
      xlab="Forcing",ylab="Onset date of wood formation (DOY)",
      cex.lab = 2, cex.axis = 2)
      m<-lm(dflat$phen~dflat$fu)
      if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkgreen")}
      mtext(paste("R2 = ",round(summary(m)$r.squared, digits = 2), sep = ""), side = 1,adj = 1, line = -2)
      mtext("D)", line = 2, adj= 0)
 
 plot(dflat$mat,dflat$phen, pch=16, bty = "l",col = "darkgray",
      xlab = "Mean Annual Temperature (C)", ylab = "Onset date of wood formation (DOY)",
      cex.lab = 2, cex.axis = 2)
     m<-lm(dflat$phen~dflat$mat)
     if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2, col = "blue")}
 
     mtext("E)", line = 2, adj= 0)
  
  latcols <-brewer.pal(5, "RdYlBu")
  
 plot(dflat$photo, dflat$phen,pch=16,col=latcols[as.factor(dflat$lat)], 
      bty="l",xlab="Photoperiod (hrs)",ylab="Onset date of wood formation (DOY)",
      cex.lab = 2, cex.axis = 2)
 m<-lm(dflat$phen~dflat$photo)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2)}
 mtext("F)", line = 2, adj= 0)
 for (l in 1:length(latz)){
   latz[l]
   dl<-as.data.frame(daylength(latz[l],-120,all_doy,8))$daylen#longitude, timezome for Seattle
   
   lines(dl, all_doy, col = alpha(latcols[l], alpha = .2), lwd = 2)
  legend("bottomright", legend = c(latz), pch = 16, col = latcols, bty = "n") 
 }
 
 
 dev.off()
 # plot(dflat$photo, dflat$mat,pch=16,col="gray", bty="l",xlab="Daylength (hrs)",ylab="Mean Annual Temperature")
 # m<-lm(dflat$mat~dflat$photo)
 # if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="gray")}
   
 # plot(dflat$lat, dflat$mat,pch=16,col="gray", bty="l",xlab="Latitude (degrees)",ylab="Mean Annual Temperature (C)")
 # m<-lm(dflat$mat~dflat$lat)
 # if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="gray")}
 # plot(dflat$lat, dflat$photo,pch=16,col="gray", bty="l",xlab="Latitude (degrees)",ylab="Photoperiod (hrs)")
 # m<-lm(dflat$mat~dflat$photo)
 # if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="gray")}

#check that temp is changing with latitude:
aggregate(dflat$mat, by=list(dflat$lat), mean)














#### Old code including chilling ###
cstar <- 110
fstaradjforchill <- 3 # how much more GDD to you need based on diff from cstar at end of daystostart

df <- data.frame(lat = numeric(), site = numeric(), year = numeric(), 
                 mat = numeric(), fu = numeric(), photo = numeric(),
                 phen = numeric(), phenrand = numeric())  

for (j in 1:sitez){
  yearly_expected_temp <- rep(basetemp, yearz)
  daily_temp <- sapply(yearly_expected_temp, function(x) c(rnorm(daysperseason, basetemp - 6, sigma),
                                                           rnorm(daysperinterseason, basetemp - 4, sigma), rnorm(daysperinterseason, basetemp - 2, sigma),
                                                           rnorm(daysperseason, basetemp, sigma)))#daily mean temperature from winter to spring (warms from mean = 0 to mean = 6)
  #chill <- daily_temp #set up chilling
  #chill[(chill)< -5] <- 0
  #chill[(chill)>5] <- 0
  force <- ifelse(daily_temp>5, 28.4/(1+exp(-.185*(daily_temp-18.4))), 0)
  fureq <- c()
  phen_dat <- c()
  for (k in 1:ncol(chill)){
    # chillsum <- sapply(1:ncol(chill), function(x) (cumsum(chill[,x])))
    forcesum <- sapply(1:ncol(force), function(x) (cumsum(force[daystostart:nrow(force),x])))
    #if (chillsum[daystostart,k]>cstar) {
    fureq[k] <- fstar
  } else {
    fureq[k] <- fstar + (cstar-chillsum[daystostart,k])*fstaradjforchill
  }
  phen_dat[k] <- min(which(forcesum[,k] > fureq[k]), na.rm = TRUE)
  meanchill <- mean(chillsum[daystostart,])#why taking mean here? mean across 30 years?
  meanfstar <- mean(fureq)
  chillmet<-length(which(chillsum[daystostart,]>cstar))/yearz
  meanforcesum<- mean(forcesum)
}#k
phen_doy<-phen_dat+daystostart#convert to doy of year
#instead of using phenology that depends on chilling and GDD, what if we just make it random?

yearly_tempall <- colMeans(daily_temp)
yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:phen_dat[x], x]))
forcesum_trunc <- sapply(1:ncol(forcesum), function(x) forcesum[phen_dat[x],x])
dl<-as.data.frame(daylength(lat,-120,phen_doy,8))$daylen#longitude, timezome for Seattle
all_doy<-seq(from =daystostart, to =(nrow(forcesum)+daystostart), by = 1)
all_dl<-as.data.frame(daylength(lat,-120,all_doy,8))$daylen#

phen_doy_rand<-rnorm(yearz,mean(phen_doy),100)#
dl_rand<-as.data.frame(daylength(lat,-120,phen_doy_rand,8))$daylen#longitude, timezome for Seattle

#plot relationships       
figname<-paste("figures/onelat",lat,yearz, "years","site",j,".pdf",sep = "_")
pdf(figname,height = 4,width = 14)
par(mfrow=c(1,4))
plot(yearly_tempall,phen_doy, pch=16, bty = "l",xlab = "Temperature", ylab = "Day of Year",ylim = range(all_doy), main = paste("site",j, sep=""))
m<-lm(phen_doy~yearly_tempall)
if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2)}

plot(forcesum_trunc, phen_doy,pch=16,col="darkgreen", bty="l",xlab="GDD",ylab="Day of Year", ylim = range(all_doy),main = paste("site",j, sep=""))
m<-lm(phen_doy~forcesum_trunc)
if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkgreen")}
#lines(forcesum, all_doy)
plot(dl, phen_doy,pch=16,col="darkblue", bty="l",xlab="Daylength (hrs)",ylab="Day of Year", ylim = range(all_doy),xlim = range(all_dl),main = paste("site",j, sep=""))
m<-lm(phen_doy~dl)
if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkblue")}
lines(all_dl, all_doy, col = alpha("darkblue", alpha = .2))
plot(dl_rand, phen_doy_rand,pch=16,col="darkblue", bty="l",xlab="Daylength (hrs)",ylab="Day of Year (Random)", ylim = range(all_doy),xlim = range(all_dl),main = paste("Random phenology site",j, sep=""))
m<-lm(phen_doy_rand~dl_rand)
if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkblue")}
lines(all_dl, all_doy, col = alpha("darkblue", alpha = .2))
dev.off()

#save all the stuff we care about   
dfadd <- data.frame(lat = lat, site=j, year = seq(1:yearz), 
                    mat =  yearly_tempall, fu=forcesum_trunc, photo = dl,
                    phen=phen_doy, phenrand = phen_doy_rand)    
df <- rbind(df, dfadd)
}#j


tail(df)
