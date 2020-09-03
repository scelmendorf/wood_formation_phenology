## Started 24 Aug 2929 ##
## By Ailene ##
## Modified from Lizzie's declining sensitivity code
## This code requires a threshold amount of forcing (aka GDD, aka fstar) to estimate phenology (e.g. wood formation onset)#
## It varies what the forcing threshold is with latitude, using a relationship estimated in Huang et al 2020 
## Forcing always starts accumulating on a certain day (Jan 1, as in Huang et al 2020) #
## After simulating phenology using only forcing, the code then looks at how well photoperiod does at predicting day of year of the phenological event
## Note: photoperiod is not used to simulate phenology doy data at all 
## We use the MAT as a starting point to simulate daily winter-spring temeprature
## We use the number of sites, years, and latitudes from the Huang et al dataset


## housekeeping

rm(list=ls()) 
options(stringsAsFactors = FALSE)

#load librarys
library(insol)
library(scales)
library(RColorBrewer)
library(dplyr)
library(colorRamps)

## Setting working directory. 
if(length(grep("ailene", getwd()))>0) { 
  setwd("C:/Users/ailene.ettinger/Documents/GitHub/wood_formation_phenology")
} 
## Read in the Huang et al dataset

d<-read.csv("data/pnas.2007058117.sd01.csv", header= TRUE)
dsub<-subset(d,select=c(site,lat,lon,year,MAT))
dsub<-distinct(dsub)#keep only distinct sites/years/mats (removing species var for now)
nyr<-rowSums(table(dsub$site,dsub$year))

## Use the mean annual temperature from Huang et al to simulate daily temperature and phenology data

daysperseason <- 60
daysperinterseason <- 30
sitez <- unique(dsub$site)

latz<-unique(dsub$lat)
lonz<-unique(dsub$lon)

#fstar <- mean(d$Force5s) #172.23 across all latitudes

#Instead of mean, let fstar vary with lat- Force5s (in their data) depends on latitude
forcelat<-aggregate(d$Force5s, by = list(d$lat), mean)
colnames(forcelat)<-c("lat","force5s")
flatmod<-lm(forcelat$force5s~forcelat$lat)

sigma <- 8

#set up colors for plotting, based on latitude
latr<-sort(unique(round(latz), digits =0))
myPalette <- colorRampPalette(brewer.pal(9, "RdYlBu")) 

latcols = myPalette(length(latr))
coltab=as.data.frame(cbind(latcols,latr))

#choose latitude to plot
lat.to.plot = "48"
col.to.plot =coltab$latcols[coltab$latr==lat.to.plot]
figname<-paste("plots/all.lats.huang.png")
png(figname,height = 800,width = 1400,res = 100)
par(mfrow=c(2,3),oma = c(1,4,1,1), mar = c(4,5,4,2))


#prepare dataframe to summarize climate, phenology across latitudes and sites
dflat <- data.frame(lat = numeric(), site = numeric(), year = numeric(), 
                    spat = numeric(), fu = numeric(), photo = numeric(),
                    phen = numeric())  
 
for (l in 1:length(latz)){
  dlat<-dsub[dsub$lat==latz[l],]
  sitelat<-unique(dlat$site)
  fstar<-coef(flatmod)[1]+coef(flatmod)[2]*latz[l]
   
  for (j in 1:length(sitelat)){
    dsite<-dlat[dlat$site == sitelat[j],]
    yearz<-dsite$year
    meantemp<- dsite$MAT
    yearly_expected_temp <- meantemp
    daily_temp <- sapply(yearly_expected_temp, function(x) c(rnorm(daysperseason, meantemp-3, sigma),
           rnorm(daysperinterseason, meantemp - 2, sigma), rnorm(daysperinterseason, meantemp, sigma),
           rnorm(daysperseason+50, meantemp+2, sigma)))#daily mean temperature from winter to spring (warms from mean = 0 (or basetemp - 6) to mean = 6 (or whatever basetemp is))

    force = ifelse(daily_temp>5, 28.4/(1+exp(-.185*(daily_temp-18.4))),0)
       
    fureq <- c()
    phen_dat <- c()
    phen_jday<-c()
       
    for (k in 1:length(yearz)){
      forcesum <- sapply(1:ncol(force), function(x) (cumsum(force[1:nrow(force),x])))
      fureq[k] <- fstar
      phen_dat[k] <- min(which(forcesum[,k] > fureq[k]))
      date<-as.character(as.Date(paste(phen_dat[k],yearz[k], sep = "-"),"%j-%Y"))
      phen_jday[k]<-JDymd(as.numeric(substr(date,1,4)), as.numeric(substr(date,6,7)),as.numeric(substr(date,9,10)))
      meanfstar <- mean(fureq)
      meanforcesum<- mean(forcesum)
      }#k
    
    yearly_tempall <- colMeans(daily_temp)#mean winter-spring temp
    yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:phen_dat[x], x]))#mean temp, ending on day of onset
    forcesum_trunc <- sapply(1:ncol(forcesum), function(x) forcesum[phen_dat[x],x])
    all_doy<-seq(from = 1, to = dim(daily_temp)[1], by = 1)
    phen_doy<-phen_dat
    dl<-as.data.frame(daylength(latz[l],lonz[l],phen_jday,8))$daylen#longitude, timezome for Seattle
    if( l==11 & sitelat[j] =="SIM"){
      plot(forcesum_trunc, phen_doy,pch=21,col = "darkgray", cex = .5,bg=col.to.plot, 
           xlim = range(forcesum), ylim = range(all_doy),bty="l",
           xlab="Forcing (Degree Days)",ylab="Onset of wood formation (DOY)",
           cex.lab = 2, cex.axis = 2)
      for (y in 1:length(yearz)){
        lines(forcesum[,y],all_doy, col = "darkgray", lwd = 1.5)}
      abline(v= fstar, lty = 2, lwd = 2,)
      points(forcesum_trunc, phen_doy,pch=21,cex = 2, col = "darkgray", bg = col.to.plot)
      
      text(fstar,max(all_doy),"FU", cex = 2)
      
      mtext("A)", line = 2, adj= 0, cex = 1.8)
    }
    
    #write out csv of daily temp and forcing by year
    colnames(daily_temp)<-colnames(forcesum)<-yearz
    rownames(daily_temp)<-rownames(forcesum)<-all_doy
    
    dtname<-paste("output/daily_temp/lat",latz[l],"site",sitelat[j],"dailytemp.csv",sep =".")
    write.csv(daily_temp,dtname)
    fname<-paste("output/forcesum/lsat",latz[l],"site",sitelat[j],"forcesum.csv",sep =".")
    write.csv(forcesum,fname)
  #save all the stuff we care about   
  dflatadd <- data.frame(lat = latz[l], site=sitelat[j], year = yearz, 
                      spat =  yearly_tempall, fu=forcesum_trunc, photo = dl,
                      phen=phen_doy)    
  dflat <- rbind(dflat, dflatadd)
  }#j
}#l
write.csv(dflat,"output/simphen.csv", row.names = FALSE)

dflat$latr<-round(dflat$lat, digits =0)

latdat<-dflat[dflat$latr == lat.to.plot,]
latcol = coltab$latcols[which(coltab$latr ==lat.to.plot)]

 plot(latdat$spat,latdat$phen, pch=21, bty = "l",cex = 2, col = "darkgray", bg = latcol,
      xlab = expression(paste("Mean Annual Temperature (",degree,"C)")), ylab = "Onset of wood formation (DOY)",
      cex.lab = 2, cex.axis = 2)
     m<-lm(latdat$phen~latdat$spat)
     if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2, col)}
     r2<-round(summary(m)$r.squared, digits = 2)
     mtext(bquote(paste('R'^'2'*' = ',.(r2))), side = 3,adj = 1, line = -3)
     mtext("B)", line = 2, adj= 0, cex = 1.8)
     
 plot(latdat$photo, latdat$phen,pch=21,cex = 2, col = "darkgray", bg=latcol,
      bty="l",xlab="Photoperiod (hrs)",ylab="Onset of wood formation (DOY)",
      cex.lab = 2, cex.axis = 2)
 m<-lm(latdat$phen~latdat$photo)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2)}
 r2<-round(summary(m)$r.squared, digits = 2)
 mtext(bquote(paste('R'^'2'*' = ',.(r2))), side = 1,adj = 1, line = -3)
 mtext("C)", line = 2, adj= 0, cex = 1.8)
 
 
 
 #now across latitudes
 coltab$latr<-as.numeric(coltab$latr)
 dflat<-left_join(dflat,coltab, by = "latr",copy = TRUE)
 
 
 plot(dflat$fu,dflat$phen,pch=21,bty = "l",cex = 2, col = "darkgray", bg=dflat$latcols, 
      xlab="Forcing (Degree Days)",ylab="Onset of wood formation (DOY)",
      cex.lab = 2, cex.axis = 2)
  mtext("D)", line = 2, adj= 0,cex = 1.8)

 plot(dflat$spat,dflat$phen, pch=21, bty = "l",cex = 2, col = "darkgray",bg  = dflat$latcols,
      xlab = expression(paste("Mean Annual Temperature (",degree,"C)")), ylab = "Onset of wood formation (DOY)",
      cex.lab = 2, cex.axis = 2)
 m<-lm(dflat$phen~dflat$spat)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2)}
 r2<-round(summary(m)$r.squared, digits = 2)
 mtext(bquote(paste('R'^'2'*' = ',.(r2))), side = 3,adj = 1, line = -3)
 
  mtext("E)", line = 2, adj= 0, cex = 1.8)
 
 plot(dflat$photo, dflat$phen,pch=21,cex = 2, col = "darkgray", bg=dflat$latcols,
      bty="l",xlab="Photoperiod (hrs)",ylab="Onset of wood formation (DOY)",
      cex.lab = 2, cex.axis = 2)
 
 for (l in 1:length(latr)){
  
   date<-as.character(as.Date(paste(all_doy,2019, sep = "-"),"%j-%Y"))
   
   all_jday<-JDymd(as.numeric(substr(date,1,4)), as.numeric(substr(date,6,7)),as.numeric(substr(date,9,10)))
   
      dl<-as.data.frame(daylength(latr[l],-120,all_jday,8))$daylen#longitude, timezome for Seattle
   
   lines(dl, all_doy, col = alpha(latcols[l], alpha = .2), lwd = 2)
 }
 m<-lm(dflat$phen~dflat$photo)
 if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2)}
 r2<-round(summary(m)$r.squared, digits = 2)
 mtext(bquote(paste('R'^'2'*' = ',.(r2))), side = 1,adj = 1, line = -3)


 mtext("F)", line = 2, adj= 0, cex = 1.8)
 legend("topleft", legend = coltab$latr[c(1,5,9,14,18,22)], pch = 21, pt.bg = coltab$latcols[c(1,5,9,14,18,22)], cex=2) 
 
 dev.off()
 
 
 #add color bar legend
 color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
   scale = (length(lut)-1)/(max-min)
   
   dev.new(width=1.75, height=5)
   plot(c(0,1), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
   axis(2, ticks, las=1)
   for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
   }	
 }
 
 color.bar(colorRampPalette(c("#D73027", "#FCA75E", "#CEEAF3", "#4575B4"))(100), min = min(latr), max = max(latr))
 
 
 