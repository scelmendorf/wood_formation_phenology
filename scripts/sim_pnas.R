# -- TODO ----

# 1 not entirely sure what the units of their ppt are, I couldn't find it
# maybe it doesn't matter? mean annual is ~908
#  so guessing mm since "globally averaged annual precipitation over 
#land is 715 mm"

#2 sort out why I get different numbers for daylength than they do. Possibly
# some sort of rounding error in the lat/longs provided in their table?
# they are close but not exactly the same though I am using the same package
#should not matter what timezone you do and indeed it doesn't but I can't match their calcs
# for (t in seq(-12, 12, by=1)){
#   cat(t)
#   cat('\n')
#   print(insol::daylength(E$lat[1], E$lon[1], E$DOY[1], tmz = t))
#   cat('\n')
# }
# E$photoperiod[1]

# dls=data.frame(doy=c(1:365), photoperiod=NA)
# for (doy in c(1:365)){
#   dls$photoperiod[doy]=insol::daylength(lat=E$lat[1], lon=E$lon[1], doy, 1)[3]
# }
# 
# dls_2=data.frame(doy=c(1:365), photoperiod=NA)
# for (doy in c(1:365)){
#   dls_2$photoperiod[doy]=insol::daylength(lat=E$lon[1], lon=E$lat[1], doy, 1)[3]
# }
# 
# E$photoperiod[1]
# dls$photoperiod[dls$doy==E$DOY[1]] #corresponds to doy ~152 even though the phenodate is 114?
# dls_2$photoperiod[dls_2$doy==E$DOY[1]] #not just bc they swapped lat long or something obvs




# 3 work on some randomizations, need to think more on how to do it smartly
# workflow would be simulate responses
# sample (without replacement - doy from the realized doys), or over sites?
# calculate x variables
# run regressions
# test significance/and or look at % var explained vs in the real data?

# -- BACKGROUND FROM THE PAPER ----
#https://www.pnas.org/content/early/2020/08/04/2007058117.short

#Chilling temperatures, which enable plants to release from the dormant state
#(51), were calculated at each site as the number of days in which the daily
#mean temperature was between ?5 and 5 °C, based on the time span commonly
#used from November 1 in the previous year to the onset date of wood formation (4)


#Forcing units 
#Forcing was computed using a sigmoid function of the average daily
#air temperature, as follows:

#Dfu=28.4/(1+exp(-.185*Tt-18.4)) if T >5, else 0

# -- SETUP ----
rm(list=ls())

## working directory should be root of this project


#where the climate data lives if not in your wd
ERA_clim_path='data/ERA_export_pnas_clim_08_20_2020.csv'

#set to true if you need to reextract clim data from earth engine
extract_EE=FALSE

#libraries
#note some other things required, mostly the functions are named but need 
# so I can sort what is what but need to have

# lubridate
# dplyr
# insol
# lme4
# MuMIn
# partR2

# Note: to install partR2 package, use the following code (this might take a while). Pls refer to https://github.com/mastoffel/partR2 for more info about partR2 package.
# install.packages("remotes")
# remotes::install_github("mastoffel/partR2", build_vignettes = FALSE, dependencies = TRUE) #Set build_vignettes = FALSE for a quick download

library (tidyverse)

#read data from the paper
E = read.csv("data/pnas.2007058117.sd01.csv")

# -- SUBSET FOR EARTH ENGINE ERA EXTRACTION ----
## Import the data 
if (extract_EE){
E_2_earth_engine=E%>%
  select(site, lat, lon)%>%
  distinct()

write.csv(E_2_earth_engine, 'data/E_2_earth_engine.csv',
          row.names=FALSE)
}

#ran the following GEE code to get the ERA daily data for their lat longs
#The E_2_earth_engine just has the latlongs

#var my_points: Table users/scelmendorf/geolocations/E_2_earth_engine

# // Import and subset ERA5 ImageCollection
# var era5_2mt = ee.ImageCollection('ECMWF/ERA5/DAILY')
# .select('mean_2m_air_temperature', 'maximum_2m_air_temperature', 'minimum_2m_air_temperature', 
#         'total_precipitation')
# .filter(ee.Filter.date('1997-01-01', '2017-01-01'));
# // NB: The dataset currently only contains data till 31 July 2019!
#   
#   
# 
# // Set name just so that the output table is easier to read.
# my_points = my_points.map(function(feature){
#   return feature.set('name', feature.id());
# });
# 
# 
# // Extract data for (Point) FeatureCollection by maping reduceRegions()
# // over the ERA5 ImageCollection. I use the first() reducer instead of the
# // mean as this makes more sense for point data (as opposed to polygons).
# var ERA5_temp_export = era5_2mt.map(function(image){
#   var temp_data = image.reduceRegions(my_points, ee.Reducer.first());
#   return(temp_data);
# }).flatten(); 
# 
# // Export FC as table to drive. 
# Export.table.toDrive({
#   collection: ERA5_temp_export,
#   description: 'ERA_export_pnas_clim_08_20_2020',
#   fileNamePrefix: 'ERA_export_pnas_clim_08_20_2020',
#   fileFormat: 'csv'
# });
# 




# -- SUBSET FOR EARTH ENGINE ERA EXTRACTION ----

#simulate responses
# sample (without replacement - doy from the realized doys)
# calculate x variables
# run regressions
# test significance

#Chilling temperatures, which enable plants to release from the dormant state
#(51), were calculated at each site as the number of days in which the daily
#mean temperature was between ?5 and 5 °C, based on the time span commonly
#used from November 1 in the previous year to the onset date of wood formation (4)


#Forcing units 
#Forcing was computed using a sigmoid function of the average daily
#air temperature, as follows:

#Dfu=28.4/(1+exp(-.185*Tt-18.4)) if T >5, else 0

# -- READ IN AND FORMAT THE CLIM DATA  ----

ERA_clim=read.csv(ERA_clim_path)

#calculate MAT off the ERA data
# use that as a bias correction for ERAs daily temp data
#ERA is a little biased high - maybe these more mountainous than grid cell avg
ERA_MAT=ERA_clim%>%
  mutate(date=lubridate::ymd(paste0(substr(system.index,1,4), '-', substr(system.index,5,6),
                                    substr(system.index,7,8))))%>%
  mutate (year=lubridate::year(date))%>%
  select(site, year, mean_2m_air_temperature)%>%
  mutate(ERA_mean_air_temp_C_unadj=mean_2m_air_temperature-270)%>%
  group_by(site, year)%>%
  summarize(ERA_MAT=mean(ERA_mean_air_temp_C_unadj), .groups='drop')%>%
  full_join(.,
            E%>%select(site, year, MAT)%>%
              distinct())%>%
  full_join(.,
            E%>%select(site, year, MAT)%>%
              distinct()%>%
              mutate(year=year-1)%>%
              rename(MAT_next=MAT))%>%
  full_join(.,
            E%>%select(site, year, MAT)%>%
              distinct()%>%
              group_by(site)%>%
              summarize(MAT_mean=mean(MAT, na.rm=TRUE), .groups='drop'))%>%
  #est true annual ppt for offset as either that year based on data
  # next year based on data
  # or site mean based on data
  mutate(MAT=case_when(
    !is.na(MAT) ~ MAT,
    !is.na(MAT_next) ~ MAT_next,
    TRUE ~ MAT_mean))%>%
  mutate(ERA_temp_offset=ERA_MAT-MAT)

  #whi is this one blanek in MAT? Actually should be ok bc it's not a year
  #that's in the data
  

#calculate annual total ppt off the ERA data
# use that as a bias correction for ERAs daily ppt data
# for some Novs, just use the next year
# for years where there is no next year data, use the mean

#note there are some oddicites in their data, e.g.
# they forgot to fill in the ppt for all trees in this site*yr
# but most have it
#View (E%>%filter(site=='5T1'&year==2002))

ERA_anPPT=ERA_clim%>%
  mutate(date=lubridate::ymd(paste0(substr(system.index,1,4), '-', substr(system.index,5,6),
                                    substr(system.index,7,8))))%>%
  mutate (year=lubridate::year(date))%>%
  select(site, year, total_precipitation)%>%
  group_by(site, year)%>%
  #group and convert to mm
  summarize(ERA_anPPT=1000*sum(total_precipitation), .groups='drop')%>%
  full_join(.,
            E%>%select(site, year, Precip)%>%
              distinct()%>%
              filter(!is.na(Precip)))%>%
  left_join(.,
            E%>%select(site, year, Precip)%>%
              distinct()%>%
              filter(!is.na(Precip))%>%
              mutate(year=year-1)%>%
              rename(Precip_next=Precip), by=c("site"="site",
                                               "year"="year"))%>%
  full_join(.,
            E%>%select(site, year, Precip)%>%
              distinct()%>%
              group_by(site)%>%
              summarize(Precip_mean=mean(Precip, na.rm=TRUE), .groups='drop'))%>%
  #est true annual ppt for offset as either that year based on data
  # next year based on data
  # or site mean based on data
  mutate(Precip=case_when(
    !is.na(Precip) ~ Precip,
    !is.na(Precip_next) ~ Precip_next,
    TRUE ~ Precip_mean))%>%
  mutate(ERA_ppt_offset_perc=ERA_anPPT/Precip)


#define function to accumulate chill from Nov to day-of-event
calc_chill=function(phen_data, clim_data){
  #df$X[duplicated(df$X)]
  df=dplyr::full_join(phen_data, clim_data, by=c('site'='site'))%>%
    filter((doy<=doy_event&year==year_of_event)|
             (doy>304&year==(year_of_event-1)))%>%
    mutate(
      chill=ifelse(-5<ERA_av_air_temp_C&ERA_av_air_temp_C<5, 1, 0))%>%
    ungroup()%>%
    arrange(Species, site, Biome, year_of_event, tree, date)%>%
    group_by(X, Species, site, Biome, year_of_event, tree)%>%
    mutate(Chill_5_5=cumsum(chill))%>%
    ungroup()%>%
    filter(doy_event==doy&year==year_of_event)%>%
    select(-chill)
  return(df)
}
#tst=calc_chill(phen_data, clim_data)
#sanity test that you didn't break something with nas or blanks
#sum (duplicated(tst$X))

#define function to accumulate ppt and forcing from Jan 1 to day-of-event
calc_ppt_force=function(phen_data, clim_data){
  df=dplyr::full_join(phen_data, clim_data, by=c('site'='site'))%>%
    filter((doy<=doy_event&!is.na(doy_event)&year_of_event==year))%>%
    mutate(
      force=ifelse(ERA_av_air_temp_C>5, 28.4/(1+exp(-.185*(ERA_av_air_temp_C-18.4))), 0))%>%
    ungroup()%>%
    arrange(Species, site, Biome, year_of_event, tree, date)%>%
    #group_by(subsite, water_year)%>%
    group_by(X, Species, site, Biome, year_of_event, tree)%>%
    #group_by(Species, site, Biome, year, tree)%>%
    mutate(Force5s=cumsum(force),
           PrecipDOY=cumsum(total_precipitation))%>%
    ungroup()%>%
    filter(doy_event==doy&year==year_of_event)
  return(df)
}
#tst=calc_ppt_force(phen_data, clim_data)


#define function to recalc all the input variables given the climate data
# and some potentially randomized phenology data
est_clim=function(phen_data, clim_data){
  #calculate chilling
  df=calc_chill(phen_data, clim_data)%>%
    #calculate forcing and total ppt
    full_join(.,  
              calc_ppt_force(phen_data, clim_data)%>%
                select(Species, site, Biome, year_of_event, tree, Force5s, PrecipDOY))%>%
    rowwise()%>%
    #sort out which scPDSI variable you want, note this part
    # will NOT work if we randomize among site yrs bc we don't have all the site/yrs for this
    mutate(date_event=as.Date(doy_event-1, origin = paste0(year, '-01-01')))%>%
    mutate(day_of_month=lubridate::day(date_event))%>%
    #mutate(month_prior_event=lubridate::month(date_event)-1)%>%
    mutate(month_prior_event=ifelse(day_of_month<=15,
                                    lubridate::month(date_event)-1, lubridate::month(date_event)))%>%
    rowwise()%>%
    mutate(photoperiod=insol::daylength(lat, lon, doy_event, 1)[3])%>%
    mutate(photoperiod=ifelse(is.finite(photoperiod), photoperiod, 24))%>%
    #sce thought they took the PDSI for the month prior but they took month of
    # you can see by looking at the actual data for 2nd half, pre month for first half
    mutate(scPDSI_p=case_when(
      month_prior_event == 1 ~ scPDSI1,
      month_prior_event == 2 ~ scPDSI2,
      month_prior_event == 3 ~ scPDSI3,
      month_prior_event == 4 ~ scPDSI4,
      month_prior_event == 5 ~ scPDSI5,
      month_prior_event == 6 ~ scPDSI6,
      TRUE ~ -9999
    ))%>%
    # #surely there is a better way to deal with the type conversions here
    mutate(scPDSI_p=ifelse(scPDSI_p==-9999, NA, scPDSI_p))#%>%
   #this is how I thought they calculated it but this is wrong
    # mutate(month_event=lubridate::month(date_event))%>%
    # mutate(scPDSI_p_2=case_when(
    #   month_event == 1 ~ scPDSI1,
    #   month_event == 2 ~ scPDSI2,
    #   month_event == 3 ~ scPDSI3,
    #   month_event == 4 ~ scPDSI4,
    #   month_event == 5 ~ scPDSI5,
    #   month_event == 6 ~ scPDSI6,
    #   TRUE ~ -9999
    # ))%>%
    # mutate(scPDSI_p_2=ifelse(scPDSI_p_2==-9999, NA, scPDSI_p_2))
  return(df)
}

ERA_clim=ERA_clim%>%
  mutate(date=lubridate::ymd(paste0(substr(system.index,1,4), '-', substr(system.index,5,6),
                                    substr(system.index,7,8))))%>%
  mutate (year=lubridate::year(date))%>%
  select(site, year, date, mean_2m_air_temperature, total_precipitation)%>%
  distinct()%>%
  rename(total_precipitation_unadj=total_precipitation)%>%
  mutate (ERA_av_air_temp_C_unadj=mean_2m_air_temperature-270)%>%
  left_join(., ERA_MAT%>%select(site, year, ERA_temp_offset), by=c('site'='site', 'year'='year'))%>% #ok
  left_join(., ERA_anPPT%>%select(site, year, ERA_ppt_offset_perc), by=c('site'='site', 'year'='year'))%>%
    mutate(ERA_av_air_temp_C=ERA_av_air_temp_C_unadj-ERA_temp_offset,
        total_precipitation=total_precipitation_unadj/ERA_ppt_offset_perc)

#if this is not broken, you should still have 365 d in this year
#View (ERA_clim%>%filter(site=='5T1'&year==2001))

# -- READ IN THE REAL DATA AND MAKE CLIMATE DATA THAT MATCHES IT  ----
# some noise here bc I used gridded climate data but mostly makes sense
# there are some errors in (I think their?) photoperiod calculations 
# in the dataset. 
phen_data=E%>%rename(year_of_event=year, doy_event=DOY)%>%
  #get rid of thier calcs (though we should at some point check these aren't
  #ridonc off what we get using ERA as the est)
select(-Chill0_5, -Chill_5_5,    
       -Chill_5_0, -Chill_10_0, -Force5s, -photoperiod, -PrecipDOY, -scPDSI_p)


clim_data=ERA_clim%>%select(site, year, date, ERA_av_air_temp_C, total_precipitation)%>%
  mutate(doy=lubridate::yday(date))

#generate our best approximation of their climate data using gridded
real_data=est_clim(phen_data, clim_data)

#sanity checks that we calculated right
#our guesstimates of chilling/forcingpretty close not exact but not super biased
plot (real_data$Force5s[order(real_data$X)]~E$Force5s[order(E$X)])
abline(0,1)

plot (real_data$Chill_5_5[order(real_data$X)]~E$Chill_5_5[order(E$X)])
abline(0,1)

#On the other hand, something is wrong w their photoperiod calculations? or mine but I did double check
plot (real_data$photoperiod[order(real_data$X)]~E$photoperiod[order(E$X)])
abline(0,1)

#they definitly have photoperiod as a function of site and doy of the event
# but I can't get the same numbers
# View (E%>%filter(site==E$site[1])%>%
#         arrange(DOY, photoperiod))

#scPDSIs now calculated the same
plot (real_data$scPDSI_p[order(real_data$X)]~E$scPDSI_p[order(E$X)])
abline(0,1)


# -- RANDOMIZE THEIR DATA IN VARIOUS WAYS AND SEE HOW THE PATTERNS LOOK ----
#first fake data - just randomize the phenodates and ignore the 
# species/site structure
out_fake=est_clim(phen_data%>%
                    mutate(doy_event=sample(phen_data$doy_event, replace=FALSE)), clim_data)
#just fakety fake data here
plot(real_data$doy_event[order(real_data$X)]~out_fake$doy_event[order(out_fake$X)])

# or randomize sites, as in, take the real phen data, strip off site, lat, lon (used for photoperd)
#and assign each site to a different site, lat and long
# does not randomize biomes of MATs or alts of sites, which is a bit of a problem
# for regenerating climate data (force, chill, MAT are wrong in this, but the
#lats/longs move so the photoperiod is faked correctly)

randomize_locs=function(df){
  site=df%>%select(site, lat, lon)%>%
    distinct()
  site$site_random=sample(site$site, replace=FALSE)
  site=left_join(site,
      site%>%rename(lat_random=lat, lon_random=lon)%>%
        select(-site_random), by=c('site_random'='site'))
  df=df%>%left_join(., site)%>%
    select(-site, -lat, -lon)%>%
    rename(site=site_random, lat=lat_random, lon=lon_random)
  return(df)
}

#prob w this is I don't have all the clim data for the matching years
#do twice just for fun, eventually make a loop/function to do these
#I think this work but should double check that it doesn't care that
# we are moving sites around between biomes this way
#or otherwise do somethign odd

##this works fine to get faked photoperiod information
#but it's suboptimal for the climate bc it needs the MAT
# to figure out the bias in the ERA clim data
# and also we don't have the droughty data for any other
# yrs at these sites
out_fake2=randomize_locs(phen_data)%>%est_clim(., clim_data)
out_fake3=randomize_locs(phen_data)%>%est_clim(., clim_data)
plot.new()
plot(real_data$doy_event[order(real_data$site, real_data$year, real_data$Species, real_data$tree)]~
       out_fake2$doy_event[order(out_fake2$site, out_fake2$year, out_fake2$Species, out_fake2$tree)])
plot(real_data$doy_event[order(real_data$site, real_data$year, real_data$Species, real_data$tree)]~
       out_fake3$doy_event[order(out_fake3$site, out_fake3$year, out_fake3$Species, out_fake3$tree)])


#functions to calculate summary statistics
photoperiod_model=function(df){
  #f1=lme4::lmer(DOY ~ photoperiod+ (1|Species)+(1|site), data = df)
  f1=lme4::lmer(doy_event ~ photoperiod + (1|Species)+(1|site), data = df)
  Var <- c(MuMIn::r.squaredGLMM(f1)[[1]], MuMIn::r.squaredGLMM(f1)[[2]]-
             MuMIn::r.squaredGLMM(f1)[[1]], 1-MuMIn::r.squaredGLMM(f1)[[2]])
  names(Var) <- c( "Photoperiod", "Random", "Unexplained")
  Var=Var*100 # Variance explained (%)
  return(list(AIC=AIC(f1), BIC=BIC(f1), Var=Var))
}

#there is actually signficantly LESS variation 
#explained by photoperiod than you would think by chance
#but this is indeed a bit weird bc I scrambled all the random effects
# sites and species and the like
photoperiod_model(real_data)

#lmer does not super love this, get some singular fits - maybe bc it's all circular? 
photoperiod_model(out_fake)

#just putting trees in random locations not their actual sites
# and calculating the daylengths from there you
# get about the same level of variation explained as you do with the real data
photoperiod_model(out_fake2)
photoperiod_model(out_fake3)

#this to me indicates their primary conclusions are not right
#Our LMM results show that the changes in the dates of wood formation onset in
#the Northern Hemisphere conifers can be modeled as a function of photoperiod
#alone, along with the random effects of site and species. The marginal R2
#and conditional R2 are 0.46 and 0.98, respectively (SI Appendix, Table S2).



mat_model=function(df){
  f1=lme4::lmer(doy_event ~ MAT+ (1|Species)+(1|site), data = df)
  #f1=lme4::lmer(DOY ~ MAT+ (1|site), data = df)
  Var <- c(MuMIn::r.squaredGLMM(f1)[[1]], MuMIn::r.squaredGLMM(f1)[[2]]-
             MuMIn::r.squaredGLMM(f1)[[1]], 1-MuMIn::r.squaredGLMM(f1)[[2]])
  names(Var) <- c("MAT", "Random", "Unexplained")
  Var=Var*100 # Variance explained (%)
  return(list(AIC=AIC(f1), BIC=BIC(f1), Var=Var))
}

#as expected the MAT doesn't show spurious correlations if you fake the data
mat_model(real_data)
mat_model(out_fake)

#the way I randomized by sites, the MAT doesn't move with the data so it wouldn't make sense to do that

force_model=function(df){
  f1=lme4::lmer(doy_event ~ Force5s+ (1|Species)+(1|site), data = df)
  Var <- c(MuMIn::r.squaredGLMM(f1)[[1]], MuMIn::r.squaredGLMM(f1)[[2]]-
             MuMIn::r.squaredGLMM(f1)[[1]], 1-MuMIn::r.squaredGLMM(f1)[[2]])
  names(Var) <- c("Force5s", "Random", "Unexplained")
  Var=Var*100 # Variance explained (%)
  return(list(AIC=AIC(f1), BIC=BIC(f1), Var=Var))
}

# in random data (maybe bc there's not so much ranefs explained partially)
# the forcing effect explains MORE of the variation
force_model(real_data)
force_model(out_fake)

chill_model=function(df){
  f1=lme4::lmer(doy_event ~ Chill_5_5+ (1|Species)+(1|site), data = df)
  Var <- c(MuMIn::r.squaredGLMM(f1)[[1]], MuMIn::r.squaredGLMM(f1)[[2]]-
             MuMIn::r.squaredGLMM(f1)[[1]], 1-MuMIn::r.squaredGLMM(f1)[[2]])
  names(Var) <- c("Chill_5_5", "Random", "Unexplained")
  Var=Var*100 # Variance explained (%)
  return(list(AIC=AIC(f1), BIC=BIC(f1), Var=Var))
}

# in random data (maybe bc there's not so much ranefs explained partially)
# the chilling effect explains MORE of the variation
chill_model(real_data)
chill_model(out_fake)


# -- RANDOMIZE THEIR DATA KEEPING SITE/SPP SUBSTRUCTURE? ----
#thinking maybe the right way to make fake data is use the real data to figure
# out the sigma_spp and sigma_sites in doys
# then just make up data from there and fit?

#to actually do this model, you want to flip it, i.e. is
#there less variation in photoperiod at date of phenology event
#than expected by chance



