############################################################################################
## R code to replicate the primary analysis, and generate results in the tables and 
## figures presented in the paper.
## Analyses were conducted with R version 3.6.2.
############################################################################################

# Clear field
rm(list=ls(all.names=T))
C:\Users\sael7303\Desktop

## Set the working directory (both Rcode.r and Data.csv are saved in this directory)
#setwd("D:\\xylem\\analysis 2020\\Rcode PNAS")
setwd("C:/Users/sael7303/Desktop/pnas_response")

## Import the data 
E = read.csv("Data.csv")

## Look at the data
dim(E)
head(E)
summary(E)

## Explanations for variables in the dataset
names(E) # Variables in the dataset
# Species: Species name
# Type: early or late successional species
# site: site name
# lat, lon, alt: latitude, longitude, and altitude of each site
# Biome: Boreal, Temperate, Mediterranean or Subtropical biome
# year: year
# tree: tree ID
# DOY: Onset date of wood formation (DOY)
# MAT: Mean annual temperature
# Chill0_5, Chill_5_5, Chill_5_0, Chill_10_0: chilling at different temperature thresholds (0, 5), (-5, 5), (-5, 0), (-10, 0), respectively
# Force5s: forcing
# photoperiod: photoperiod (hrs)
# scPDSI1, scPDSI2, scPDSI3, scPDSI4, scPDSI5, scPDSI6: scPDSI (self-calibrating Palmer Drought Severity Index) in Jan, Feb, ..., June
# scPDSI_p: scPDSI in the month prior to the onset date of wood formation (DOY)
# SPEI1m1, SPEI1m2, SPEI1m3, SPEI1m4, SPEI1m5, SPEI1m6: 1-month SPEI (Standardized Precipitation Evapotranspiration Index) in Jan, Feb, ..., June
# SPEI1m_p: 1-month SPEI in the month prior to the onset date of wood formation (DOY)
# SPEI3m1, SPEI3m2, SPEI3m3, SPEI3m4, SPEI3m5, SPEI3m6: 3-month SPEI in Jan, Feb, ..., June
# SPEI3m_p: 3-month SPEI in the month prior to the onset of wood formation
# SPEI6m1, SPEI6m2, SPEI6m3, SPEI6m4, SPEI6m5, SPEI6m6: 6-month SPEI in Jan, Feb, ..., June
# SPEI6m_p: 6-month SPEI in the month prior to the onset date of wood formation (DOY)
# Precip: Annual total precipitation
# PrecipDOY: precipitation from Jan, 1 to the onset date of wood formation (DOY)

## Reorder factor levels of Biome
E$Biome = factor(E$Biome, levels = c("Boreal","Temperate", "Mediterranean", "Subtropical"))
table(E$Biome)

## Load required packages. Pls make sure all packages are installed. 
library(car)
library(lme4)
library(lmerTest)
library(ggplot2)
library(mgcv)
library(MuMIn)
library(remotes)
# Note: to install partR2 package, use the following code (this might take a while). Pls refer to https://github.com/mastoffel/partR2 for more info about partR2 package.
# install.packages("remotes")
# remotes::install_github("mastoffel/partR2", build_vignettes = FALSE, dependencies = TRUE) #Set build_vignettes = FALSE for a quick download
library(partR2) 

############################################################################################
## Code for Tables
############################################################################################

############################################################################################
## Table 1: Statistics of Model 1, including photoperiod, forcing, chilling, and scPDSI, and 
## of Model 2, including MAT, forcing, chilling, and scPDSI. Model 2 was applied to all data 
## and to each group (boreal and temperate biomes; early and late successional species)

## Note: Compared with scPDSI in other months, scPDSI_p (scPDSI in the month prior to the 
# onset date of wood formation (DOY)) had the highest R2, and thus was included in the model. 
# See the code at the very end for the comparison of different scPDSIs. The results of the 
# SPEI were similar to those for the scPDSI and are reported in the Table S4.

# Model 1: DOY ~ Photoperiod + forcing + chilling + scPDSI
f1=lmer(DOY ~ photoperiod + Force5s + Chill_5_5 + scPDSI_p + (1|Species)+(1|site), data = E)# Fitting a mixed-effects model with lmer
vif(f1) # Check multicollinearity
summary(f1) # summary of the model
r.squaredGLMM(f1)# # explained variance by each part
# R2m (marginal): variance explained by the fixed effects
# R2c (conditional): variance explained by the entire model, including both fixed and random effects
c(AIC(f1), BIC(f1)) # AIC, BIC
hist(residuals(f1)) # check normality of residuals

# Model 2: DOY ~ MAT + forcing + chilling + scPDSI
# Model 2 for all data
f2=lmer(DOY ~ MAT + Force5s + Chill_5_5 + scPDSI_p + (1|Species)+(1|site), data = E)
vif(f2)
summary(f2)
r.squaredGLMM(f2)
c(AIC(f2), BIC(f2))
qqnorm(residuals(f2)) # check normality of residuals
hist(residuals(f2))

# Model 2 for Boreal biome
f2b=lmer(DOY ~ MAT + Force5s + Chill_5_5 + scPDSI_p +(1|site)+(1|Species), data = E[E$Biome=="Boreal",])
vif(f2b)
summary(f2b)
c(AIC(f2b), BIC(f2b))
r.squaredGLMM(f2b)
qqnorm(residuals(f2b))

# Model 2 for Temperate biome
f2t=lmer(DOY ~ MAT + Force5s + Chill_5_5 + scPDSI_p +(1|site)+(1|Species), data = E[E$Biome=="Temperate",])
vif(f2t)
summary(f2t)
c(AIC(f2t), BIC(f2t))
r.squaredGLMM(f2t)
qqnorm(residuals(f2t))

# Model 2 for Early successional species
f2e=lmer(DOY ~MAT + Force5s + Chill_5_5 + scPDSI_p + (1|Species)+(1|site), data = E[E$Type=="Early",])
vif(f2e)
summary(f2e)
c(AIC(f2e), BIC(f2e))
r.squaredGLMM(f2e)
qqnorm(residuals(f2e))

# Model 2 for Late successional species
f2l=lmer(DOY ~MAT + Force5s + Chill_5_5 + scPDSI_p + (1|Species)+(1|site), data = E[E$Type=="Late",])
vif(f2l)
summary(f2l)
c(AIC(f2l), BIC(f2l))
r.squaredGLMM(f2l)
qqnorm(residuals(f2l))

############################################################################################
## Table S2: Statistics of Models: DOY ~ Photoperiod, and DOY ~ MAT
# Statistics of model: DOY ~ Photoperiod
f=lmer(DOY ~ photoperiod+ (1|Species)+(1|site), data = E)
summary(f)
r.squaredGLMM(f)
c(AIC(f), BIC(f))
Var <- c(r.squaredGLMM(f)[[1]], r.squaredGLMM(f)[[2]]-r.squaredGLMM(f)[[1]], 1-r.squaredGLMM(f)[[2]])
names(Var) <- c( "Photoperiod", "Random", "Unexplained")
Var*100 # Variance explained (%)

# Statistics of model: DOY ~ MAT
f=lmer(DOY ~ MAT+ (1|Species)+(1|site), data = E)
summary(f)
r.squaredGLMM(f) 
c(AIC(f), BIC(f)) 
Var <- c(r.squaredGLMM(f)[[1]], r.squaredGLMM(f)[[2]]-r.squaredGLMM(f)[[1]], 1-r.squaredGLMM(f)[[2]])
names(Var) <- c( "MAT", "Random", "Unexplained")
Var*100 # Variance explained (%)

############################################################################################
## Table S3: Statistics of Model: DOY ~ MAT + Forcing + Chilling + Precipitation
# Statistics of model using annual total precipitation
f=lmer(DOY ~ MAT + Force5s + Chill_5_5 +Precip+ (1|Species)+(1|site), data = E[!is.na(E$Precip),])
summary(f)
r.squaredGLMM(f)
c(AIC(f), BIC(f))
hist(residuals(f)) 
qqnorm(residuals(f)) 

# Statistics of model using precipitation (Jan, 1st to DOY)
f=lmer(DOY ~ MAT + Force5s + Chill_5_5 +PrecipDOY+ (1|Species)+(1|site), data = E[!is.na(E$PrecipDOY),])
summary(f)
r.squaredGLMM(f)
c(AIC(f), BIC(f))
hist(residuals(f)) 
qqnorm(residuals(f)) 

# Variance partition of the studied variables
# Function var_decomp: Given lmer model, variable names, and dataset, var_decomp partitions the variances
# explained by each fixed variable, random variable, and unexplained component.
var_decomp <- function(f, var, dat){
  R2<-partR2(f, partvars = var, R2_type = "marginal",
             max_level = 1, nboot = 10, CI = 0.95, data = dat) # compute partial R2
  PR2_fixed = R2$R2[c(2,3,4,5),2] # extract the partial R2 for four fixed variables
  PR2_fixed[PR2_fixed<0]=0 # if the contribution of a variable is too little, small negative values may occur and were set to 0.
  PR2_fixed = PR2_fixed/sum(PR2_fixed)*r.squaredGLMM(f)[[1]] # rescale partial R2 of fixed variables to sum to marginal R2
  PR2_all <- c(PR2_fixed, r.squaredGLMM(f)[[2]]-r.squaredGLMM(f)[[1]], 1-r.squaredGLMM(f)[[2]]) # combine partial R2 for four fixed variables, random effects, and unexplained variance
  return(PR2_all*100) # return partitioned R2 values (%)
}
var = c("MAT","Force5s","Chill_5_5", "PrecipDOY")
PR2_all = var_decomp(f, var, E[!is.na(E$PrecipDOY),]) # call function var_decomp to compute partitioned variances
names(PR2_all) <- c( var, "Random", "Unexplained")
PR2_all # Variances explained by each fixed variable, random variable, and unexplained component

############################################################################################
## Table S4: Comparison among models based on different SPEI temporal scales: 1-month, 3-month and 6-month. 
## Statistics of Model: DOY ~ MAT + Forcing + Chilling + SPEI
# Model using 1-month SPEI
f1m=lmer(DOY ~ MAT + Force5s + Chill_5_5 +SPEI1m_p+ (1|Species)+(1|site), data = E)
summary(f1m)
r.squaredGLMM(f1m)
c(AIC(f1m), BIC(f1m))
hist(residuals(f1m)) 
# Model using 3-month SPEI
f3m=lmer(DOY ~ MAT + Force5s + Chill_5_5 +SPEI3m_p+ (1|Species)+(1|site), data = E)
summary(f3m)
r.squaredGLMM(f3m)
c(AIC(f3m), BIC(f3m))
hist(residuals(f3m)) 
# Model using 6-month SPEI
f6m=lmer(DOY ~ MAT + Force5s + Chill_5_5 +SPEI6m_p+ (1|Species)+(1|site), data = E)
summary(f6m)
r.squaredGLMM(f6m)
c(AIC(f6m), BIC(f6m))
hist(residuals(f6m)) 

############################################################################################
## Table S5: Model results of Mediterranean biome: DOY ~ MAT + Forcing + Chilling + scPDSI
f2m=lmer(DOY ~ MAT + Force5s + Chill_5_5 + scPDSI_p +(1|site)+(1|Species), data = E[E$Biome=="Mediterranean",])
summary(f2m)
r.squaredGLMM(f2m)
c(AIC(f2m), BIC(f2m))
# Variance partition of the studied variables
var = c("MAT","Force5s","Chill_5_5", "scPDSI_p")
PR2_all = var_decomp(f2m, var, E[E$Biome=="Mediterranean",])
names(PR2_all) <- c( var, "Random", "Unexplained")
PR2_all# Variances explained by each fixed variable, random variable, and unexplained component

############################################################################################
## Table S6: Statistics of Model 2 at different chilling thresholds.
# Chill0_5: chilling between 0 and 5 degree
f=lmer(DOY ~ MAT + Force5s + Chill0_5 + scPDSI_p + (1|Species)+(1|site), data = E)
summary(f)
r.squaredGLMM(f)
# Variance partition of the studied variables
var = c("MAT","Force5s","Chill0_5", "scPDSI_p")
PR2_all = var_decomp(f, var, E)
names(PR2_all) <- c( var, "Random", "Unexplained")
PR2_all # Variances explained by each fixed variable, random variable, and unexplained component

# Chill_5_5: chilling between -5 and 5 degree
f=lmer(DOY ~ MAT + Force5s + Chill_5_5 + scPDSI_p + (1|Species)+(1|site), data = E)
summary(f)
r.squaredGLMM(f)
# Variance partition of the studied variables
var = c("MAT","Force5s","Chill_5_5", "scPDSI_p")
PR2_all = var_decomp(f, var, E)
names(PR2_all) <- c( var, "Random", "Unexplained")
PR2_all # Variances explained by each fixed variable, random variable, and unexplained component

# Chill_5_0: chilling between -5 and 0 degree
f=lmer(DOY ~ MAT + Chill_5_0 + Force5s + scPDSI_p + (1|Species)+(1|site), data = E)
summary(f)
r.squaredGLMM(f)
# Variance partition of the studied variables
var = c("MAT","Force5s","Chill_5_0", "scPDSI_p")
PR2_all = var_decomp(f, var, E)
names(PR2_all) <- c( var, "Random", "Unexplained")
PR2_all # Variances explained by each fixed variable, random variable, and unexplained component

# Chill_10_0: chilling between -10 and 0 degree
f=lmer(DOY ~ MAT + Chill_10_0 + Force5s+ scPDSI_p  + (1|Species)+(1|site), data = E)
summary(f)
r.squaredGLMM(f)
# Variance partition of the studied variables
var = c("MAT","Force5s","Chill_10_0","scPDSI_p")
PR2_all = var_decomp(f, var, E)
names(PR2_all) <- c( var, "Random", "Unexplained")
PR2_all # Variances explained by each fixed variable, random variable, and unexplained component


############################################################################################
## Code for Figures
############################################################################################
## create a folder in current working directory to save plots
dir.create("plots")

############################################################################################
## Fig1: onset date (DOY) vs latitude
tiff(file=paste("plots/Fig1_DOY_latitude.tiff",sep=""), width = 8, height = 4, units = 'in', res = 300, compression = 'lzw') #save the plot in tiff fromat
f = gam(DOY~s(lat, k=4),data=E) #gam model
summary(f)
xbreaks = seq(20,70,10)
xlabels <- unlist(lapply(xbreaks, function(x) paste(x,"\u00B0N"))) # format latitude axis labels
ggplot() +
  geom_point(aes(lat, DOY, color=Species, shape = Biome), data = E, size = 1.5)+ theme_bw()+
  xlab("Latitude") +  ylab("Onset date of wood formation (DOY)")+
  geom_smooth(aes(lat, DOY), data=E, method = "gam",formula= y ~ s(x, k=4), alpha = 0.5, lwd = 2)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 15))+
  annotate("text", x = 50, y = 30,label = paste("EDF= 2.8, P<0.001"))+ #EDF and p-values are calculated with gam model
  scale_x_continuous(limits = c(20, 70), breaks = xbreaks, labels = xlabels)+
  guides(shape =guide_legend(ncol=2), color = guide_legend(ncol=3))
dev.off()

############################################################################################
## Fig 2a: onset date (DOY) vs MAT vs Forcing
# step 1: gam model and prepare fitted values
model = gam(DOY ~ s(MAT) + s(Force5s), data=E)
summary(model)
x = seq(-5, 25, 2) # range for MAT
y = seq(100,800,30) # range for Forcing
f <- function(x, y) {
  r = predict.gam(model, newdata=data.frame(MAT=x, Force5s=y))}
z = outer(x, y, f) # fitted DOY
nrz = nrow(z)
ncz = ncol(z)
jet.colors = colorRampPalette( c("blue","green", "yellow", "red") )# Define the range of colors
nbcol = 20 # The number of colors from this palette
color=jet.colors(nbcol)
zfacet = (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4 # Compute the z-value at the center of each facet
facetcol = cut(zfacet, nbcol) # Convert facet z-values into color indices
# step2: generate 3D plot
tiff(file=paste("plots/Fig2a_DOY_MAT_Forcing.tiff",sep=""), width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
layout(matrix(1:2, nrow=1, ncol=2), widths=c(4,1), heights=1)
par(bg = "white", mar=c(0,2,0,1))
persp(x,y,z,theta = 35, phi = 20, expand = 0.7, col = color[facetcol],
      xlab="Temperature (\u00B0C)", ylab="Forcing (FU)", zlab="Onset date of wood formation (DOY)",
      r = 10,ticktype = "detailed")
labs = levels(facetcol)
tmp = cbind(lower = as.numeric(sub("\\((.+),.*", "\\1", labs)),
            upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", labs)))
par(mar=c(7,0,7,5))
image(x=1, y=rowMeans(tmp), matrix(rowMeans(tmp), nrow=1, ncol=nbcol), col=color, axes=FALSE, xlab="", ylab="")
axis(4)
dev.off()

############################################################################################
## Fig 2b onset date (DOY) vs MAT vs chilling
# step 1: gam model and prepare fitted values
model = gam(DOY ~ s(MAT) + s(Chill_5_5), data=E)
x = seq(-5, 25, 2)
y = seq(0,160,10)
f = function(x, y) {
  r = predict.gam(model, newdata=data.frame(MAT=x, Chill_5_5=y))}
z = outer(x, y, f)
nrz = nrow(z)
ncz = ncol(z)
jet.colors = colorRampPalette( c("blue","green", "yellow", "red") )
nbcol = 20
color=jet.colors(nbcol)
zfacet = (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4
facetcol = cut(zfacet, nbcol)
# step2: generate 3D plot
tiff(file=paste("plots/Fig2b_DOY_MAT_Chilling.tiff",sep=""), width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
layout(matrix(1:2, nrow=1, ncol=2), widths=c(4,1), heights=1)
par(bg = "white", mar=c(0,2,0,1))
persp(x,y,z,theta = 35, phi = 20, expand = 0.7, col = color[facetcol],
      xlab="Temperature (\u00B0C)", ylab="Chilling (days)", zlab="Onset date of wood formation (DOY)",
      r = 10,ticktype = "detailed")
labs = levels(facetcol)
tmp = cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
            upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
par(mar=c(7,0,7,5))
image(x=1, y=rowMeans(tmp), matrix(rowMeans(tmp), nrow=1, ncol=nbcol), col=color, axes=FALSE, xlab="", ylab="")
axis(4)
dev.off()

############################################################################################
## Fig2c: onset date (DOY) vs forcing vs chilling
# step 1: gam model and prepare fitted values
model = gam(DOY ~ s(Force5s) + s(Chill_5_5), data=E)
x = seq(100,800,30)
y = seq(0,160,10)
f = function(x, y) {
  r = predict.gam(model, newdata=data.frame(Force5s=x, Chill_5_5=y))}
z = outer(x, y, f)
nrz = nrow(z)
ncz = ncol(z)
jet.colors = colorRampPalette( c("blue","green", "yellow", "red") )
nbcol = 20
color=jet.colors(nbcol)
zfacet = (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4
facetcol = cut(zfacet, nbcol)
# step2: generate 3D plot
tiff(file=paste("plots/Fig2c_DOY_Forcing_Chilling.tiff",sep=""), width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
layout(matrix(1:2, nrow=1, ncol=2), widths=c(4,1), heights=1)
par(bg = "white", mar=c(0,2,0,1))
persp(x,y,z,theta = 30, phi = 30, expand = 0.7, col = color[facetcol],
      xlab="Forcing (FU)", ylab="Chilling (days)", zlab="Onset date of wood formation (DOY)",
      r = 10,ticktype = "detailed")
labs = levels(facetcol)
tmp = cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
            upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
par(mar=c(7,0,7,5))
image(x=1, y=rowMeans(tmp), matrix(rowMeans(tmp), nrow=1, ncol=nbcol), col=color, axes=FALSE, xlab="", ylab="")
axis(4)
dev.off()

############################################################################################
## Fig 3: Compute the partitioned R2 for model 1 and model 2 in Table 1. Then barplot of partitioned R2 
## in Fig. 3 was generated in sigmaplot
# Function var_decomp is repeated here for convenience
# Function var_decomp: Given lmer model, variable names, and dataset, var_decomp partitions the variances
# explained by each fixed variable, random variable, and unexplained component
var_decomp <- function(f, var, dat){
  R2<-partR2(f, partvars = var, R2_type = "marginal",
             max_level = 1, nboot = 10, CI = 0.95, data = dat) # compute partial R2
  PR2_fixed = R2$R2[c(2,3,4,5),2] # extract the partial R2 for four fixed variables
  PR2_fixed[PR2_fixed<0]=0 # if the contribution of a variable is too little, small negative values may occur and were set to 0.
  PR2_fixed = PR2_fixed/sum(PR2_fixed)*r.squaredGLMM(f)[[1]] # rescale partial R2 of fixed variables to sum to marginal R2
  PR2_all <- c(PR2_fixed, r.squaredGLMM(f)[[2]]-r.squaredGLMM(f)[[1]], 1-r.squaredGLMM(f)[[2]]) # combine partial R2 for four fixed variables, random effects, and unexplained variance
  return(PR2_all*100) # return partitioned R2 values (%)
}

## Variance partition of the studied variables in model 1 and model 2
# Model 1: DOY ~ Photoperiod + forcing + chilling + scPDSI
f1=lmer(DOY ~ photoperiod + Force5s + Chill_5_5 + scPDSI_p + (1|Species)+(1|site), data = E)
var = c("photoperiod","Force5s","Chill_5_5", "scPDSI_p")
PR2_all = var_decomp(f1, var, E)
names(PR2_all) <- c( var, "Random", "Unexplained")
PR2_all # Variances explained by each fixed variable, random variable, and unexplained component

# Model 2: DOY ~ MAT + forcing + chilling + scPDSI
# Model 2 for all data
f2=lmer(DOY ~ MAT + Force5s + Chill_5_5 + scPDSI_p + (1|Species)+(1|site), data = E)
var = c("MAT","Force5s","Chill_5_5","scPDSI_p")
PR2_all = var_decomp(f2, var, E)
names(PR2_all) <- c( var, "Random", "Unexplained")
PR2_all # Variances explained by each fixed variable, random variable, and unexplained component

# Model 2 for Boreal biome
f2b=lmer(DOY ~ MAT + Force5s + Chill_5_5 + scPDSI_p +(1|site)+(1|Species), data = E[E$Biome=="Boreal",])
var = c("MAT","Force5s","Chill_5_5","scPDSI_p")
PR2_all = var_decomp(f2b, var, E[E$Biome=="Boreal",])
names(PR2_all) <- c( var, "Random", "Unexplained")
PR2_all # Variances explained by each fixed variable, random variable, and unexplained component

# Model 2 for Temperate biome
f2t=lmer(DOY ~ MAT + Force5s + Chill_5_5 + scPDSI_p +(1|site)+(1|Species), data = E[E$Biome=="Temperate",])
var = c("MAT","Force5s","Chill_5_5","scPDSI_p")
PR2_all = var_decomp(f2t, var, E[E$Biome=="Temperate",])
names(PR2_all) <- c( var, "Random", "Unexplained")
PR2_all # Variances explained by each fixed variable, random variable, and unexplained component

# Model 2 for Early successional species
f2e=lmer(DOY ~MAT + Force5s + Chill_5_5 + scPDSI_p + (1|Species)+(1|site), data = E[E$Type=="Early",])
var = c("MAT","Force5s","Chill_5_5","scPDSI_p")
PR2_all = var_decomp(f2e, var, E[E$Type=="Early",])
names(PR2_all) <- c( var, "Random", "Unexplained")
PR2_all # Variances explained by each fixed variable, random variable, and unexplained component

# Model 2 for Late successional species
f2l=lmer(DOY ~MAT + Force5s + Chill_5_5 + scPDSI_p + (1|Species)+(1|site), data = E[E$Type=="Late",])
var = c("MAT","Force5s","Chill_5_5","scPDSI_p")
PR2_all = var_decomp(f2l, var, E[E$Type=="Late",])
names(PR2_all) <- c( var, "Random", "Unexplained")
PR2_all # Variances explained by each fixed variable, random variable, and unexplained component

############################################################################################
## Fig S2: MAT vs photoperiod
summary(lm(MAT~photoperiod, E)) # linear model
tiff(file=paste("plots/FigS2_MAT_photoperiod.tiff",sep=""), width = 8, height = 4, units = 'in', res = 300, compression = 'lzw')
ggplot() +
  geom_point(aes(photoperiod, MAT, color=Species, shape = Biome), data=E, size = 1.5)+ theme_bw()+
  xlab("Photoperiod (hrs)") +  ylab(expression("Mean annual temperature ("*degree~C *")"))+
  geom_smooth(aes(photoperiod, MAT), data=E, method = "gam",formula= y ~ x, alpha = 0.5, lwd = 2)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 15))+
  annotate("text", x = 19, y = 20, label = "paste(italic(R) ^ 2, \" = 0.617, P<0.001\")", parse = TRUE)+ # R2 and p-values are from lm
  scale_x_continuous(limits = c(10, 24), breaks = seq(10,24,2))+
  scale_y_continuous(limits = c(-3, 23), breaks = seq(0,20,5))+
  guides(shape =guide_legend(ncol=2), color = guide_legend(ncol=3))#column of legend
dev.off()

############################################################################################
## Fig S3: MAT vs latitude
tiff(file=paste("plots/FigS3_MAT_latitude.tiff",sep=""), width = 8, height = 4, units = 'in', res = 300, compression = 'lzw')
f = gam(MAT~s(lat, k=4),data=E) #gam model
summary(f)
xbreaks = seq(20,70,10)
xlabels <- unlist(lapply(xbreaks, function(x) paste(x,"\u00B0N")))
ggplot() +
  geom_point(aes(lat, MAT, color=Species, shape = Biome), data=E, size = 1.5)+ theme_bw()+
  xlab("Latitude") +  ylab(expression("Mean annual temperature ("*degree~C *")"))+
  geom_smooth(aes(lat, MAT), data=E, method = "gam",formula= y ~ s(x, k=4), alpha = 0.5, lwd = 2)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 15))+
  annotate("text", x = 50, y = 20,label = paste("EDF= 2.9, P<0.001"))+
  scale_x_continuous(limits = c(20, 70), breaks = xbreaks, labels = xlabels)+
  guides(shape =guide_legend(ncol=2), color = guide_legend(ncol=3))#column of legend
dev.off()

############################################################################################
## Fig S4: Onset date vs MAT
f = gam(DOY~s(MAT,k=4),data=E) #gam model
summary(f)
tiff(file=paste("plots/FigS4_DOY_MAT.tiff",sep=""), width = 8, height = 4, units = 'in', res = 300, compression = 'lzw')
ggplot() +
  geom_point(aes(MAT, DOY, shape = Biome, color=Species), data=E, size = 1.5)+ theme_bw()+
  xlab(expression("Mean annual temperature ("*degree~C *")")) +  ylab("Onset date of wood formation (DOY)")+
  geom_smooth(aes(MAT, DOY), data=E, method = "gam",formula= y ~ s(x,k=4), alpha = 0.5, lwd = 2)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 15))+
  annotate("text", x = 5, y = 30,label = paste("EDF= 2.9, P<0.001"))+
  scale_x_continuous(limits = c(-2, 23), breaks = seq(-10,25,5))+
  guides(shape =guide_legend(ncol=2), color = guide_legend(ncol=3))#column of legend
dev.off()

############################################################################################
## FigS5a: photoperiod vs latitude
tiff(file=paste("plots/FigS5a_photoperiod_latitude.tiff",sep=""), width = 8, height = 4, units = 'in', res = 300, compression = 'lzw')
f = gam(photoperiod~s(lat, k=4),data=E) #gam model
summary(f)
xbreaks = seq(20,70,10)
xlabels <- unlist(lapply(xbreaks, function(x) paste(x,"\u00B0N")))
ggplot() +
  geom_point(aes(lat, photoperiod, color=Species, shape = Biome), data=E, size = 1.5)+ theme_bw()+
  xlab("Latitude") +  ylab("Photoperiod (hrs)")+
  geom_smooth(aes(lat, photoperiod), data=E, method = "gam",formula= y ~ s(x, k =4), alpha = 0.5, lwd = 2)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 15))+
  annotate("text", x = 55, y = 11,label = paste("EDF= 2.9, P<0.001"))+
  scale_x_continuous(limits = c(20, 70), breaks = xbreaks, labels = xlabels)+
  scale_y_continuous(limits = c(10, 24), breaks = seq(10,24,2))+
  guides(shape =guide_legend(ncol=2), color = guide_legend(ncol=3))#column of legend
dev.off()

############################################################################################
## FigS5b: onset date vs photoperiod
summary(lm(DOY~photoperiod, E))
tiff(file=paste("plots/FigS5b_DOY_photoperiod.tiff",sep=""), width = 8, height = 4, units = 'in', res = 300, compression = 'lzw')
ggplot() +
  geom_point(aes(photoperiod, DOY, color=Species, shape = Biome), data=E, size = 1.5)+ theme_bw()+
  xlab("Photoperiod (hrs)") +  ylab("Onset date of wood formation (DOY)")+
  geom_smooth(aes(photoperiod, DOY), data=E, method = "gam",formula= y ~ x, alpha = 0.5, lwd = 2)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 15))+
  annotate("text", x = 15, y = 11,label = "paste(italic(R) ^ 2, \" = 0.677, P<0.001\")", parse = TRUE)+
  scale_x_continuous(limits = c(10, 24), breaks = seq(10,24,2))+
  scale_y_continuous(limits = c(0, 220), breaks = seq(0,220,50))+
  guides(shape =guide_legend(ncol=2), color = guide_legend(ncol=3))#column of legend
dev.off()

############################################################################################
## Fig S6: onset date (DOY) vs MAT vs scPDSI
# step 1: gam model and prepare fitted values
model = gam(DOY ~ s(MAT) + s(scPDSI_p), data=E)
summary(model)
x = seq(-5, 25, 2)
y = seq(-4.5, 4.5, 0.5)
f = function(x, y) {
  r = predict.gam(model, newdata=data.frame(MAT=x, scPDSI_p=y))}
z = outer(x, y, f)
nrz = nrow(z)
ncz = ncol(z)
jet.colors = colorRampPalette( c("blue","green", "yellow", "red") )
nbcol = 20
color=jet.colors(nbcol)
zfacet = (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4
facetcol = cut(zfacet, nbcol)
# step2: generate 3D plot
tiff(file=paste("plots/FigS6_DOY_MAT_scPDSI.tiff",sep=""), width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
layout(matrix(1:2, nrow=1, ncol=2), widths=c(4,1), heights=1)
par(bg = "white", mar=c(0,2,0,1))
persp(x,y,z,theta = 35, phi = 20, expand = 0.7, col = color[facetcol],
      xlab="Temperature (\u00B0C)", ylab="scPDSI", zlab="Onset date of wood formation (DOY)",
      r = 10,ticktype = "detailed")
labs = levels(facetcol)
tmp = cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
             upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
par(mar=c(7,0,7,5))
image(x=1, y=rowMeans(tmp), matrix(rowMeans(tmp), nrow=1, ncol=nbcol), col=color, axes=FALSE, xlab="", ylab="")
axis(4)
dev.off()

############################################################################################
## Fig S7a: onset date (DOY) vs Photoperiod vs Forcing
# step 1: gam model and prepare fitted values
model = gam(DOY ~ s(photoperiod) + s(Force5s), data=E)
summary(model)
x = seq(10, 24, 1)
y = seq(100,800,30)
f = function(x, y) {
  r = predict.gam(model, newdata=data.frame(photoperiod=x, Force5s=y))}
z = outer(x, y, f)
nrz = nrow(z)
ncz = ncol(z)
jet.colors = colorRampPalette( c("blue","green", "yellow", "red") )
nbcol = 20
color=jet.colors(nbcol)
zfacet = (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4
facetcol = cut(zfacet, nbcol)
# step2: generate 3D plot
tiff(file=paste("plots/FigS7a_DOY_photoperiod_Forcing.tiff",sep=""), width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
layout(matrix(1:2, nrow=1, ncol=2), widths=c(4,1), heights=1)
par(bg = "white", mar=c(0,2,0,1))
persp(x,y,z,theta = -60, phi = 20, expand = 0.7, col = color[facetcol],
             xlab="Photoperiod (hrs)", ylab="Forcing (FU)", zlab="Onset date of wood formation (DOY)",
             r = 10,ticktype = "detailed")
labs = levels(facetcol)
tmp = cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
             upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
par(mar=c(7,0,7,5))
image(x=1, y=rowMeans(tmp), matrix(rowMeans(tmp), nrow=1, ncol=nbcol), col=color, axes=FALSE, xlab="", ylab="")
axis(4)
dev.off()

############################################################################################
## Fig S7b: onset date (DOY) vs Photoperiod vs chilling
# step 1: gam model and prepare fitted values
model = gam(DOY ~ s(photoperiod) + s(Chill_5_5), data=E)
x = seq(10, 24, 1)
y = seq(0,160,10)
f = function(x, y) {
  r = predict.gam(model, newdata=data.frame(photoperiod=x, Chill_5_5=y))}
z = outer(x, y, f)
nrz = nrow(z)
ncz = ncol(z)
jet.colors = colorRampPalette( c("blue","green", "yellow", "red") )
nbcol = 20
color=jet.colors(nbcol)
zfacet = (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4
facetcol = cut(zfacet, nbcol)
# step2: generate 3D plot
tiff(file=paste("plots/Figs7b_DOY_photoperiod_Chilling.tiff",sep=""), width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
layout(matrix(1:2, nrow=1, ncol=2), widths=c(4,1), heights=1)
par(bg = "white", mar=c(0,2,0,1))
persp(x,y,z,theta = -60, phi = 20, expand = 0.7, col = color[facetcol],
      xlab="Photoperiod (hrs)", ylab="Chilling (days)", zlab="Onset date of wood formation (DOY)",
      r = 10,ticktype = "detailed")
labs = levels(facetcol)
tmp = cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
             upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
par(mar=c(7,0,7,5))
image(x=1, y=rowMeans(tmp), matrix(rowMeans(tmp), nrow=1, ncol=nbcol), col=color, axes=FALSE, xlab="", ylab="")
axis(4)
dev.off()

############################################################################################
## Code for preliminary analysis
############################################################################################
## Compared with scPDSI in other months, scPDSI_p (scPDSI in the month prior to the 
# onset date of wood formation (DOY)) had the highest R2, and thus was included in the model.

# scPDSI in the month prior to the onset date of wood formation (DOY)
f=lmer(DOY ~ scPDSI_p + (1|Species)+(1|site), data = E) ## R2 is the highest
summary(f)   
r.squaredGLMM(f)
# scPDSI in Jan
f=lmer(DOY ~ scPDSI1 + (1|Species)+(1|site), data = E) 
summary(f)   
r.squaredGLMM(f)
# scPDSI in Feb
f=lmer(DOY ~ scPDSI2 + (1|Species)+(1|site), data = E) 
summary(f)   
r.squaredGLMM(f)
# scPDSI in Mar
f=lmer(DOY ~ scPDSI3 + (1|Species)+(1|site), data = E) 
summary(f)   
r.squaredGLMM(f)
# scPDSI in Apr
f=lmer(DOY ~ scPDSI4 + (1|Species)+(1|site), data = E) 
summary(f)   
r.squaredGLMM(f)
# scPDSI in May
f=lmer(DOY ~ scPDSI5 + (1|Species)+(1|site), data = E) 
summary(f)   
r.squaredGLMM(f)
# scPDSI in Jun
f=lmer(DOY ~ scPDSI6 + (1|Species)+(1|site), data = E) 
summary(f)   
r.squaredGLMM(f)

