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
#   dls$photoperiod[doy]=insol::daylength(lat=48.35, lon=E$lon[1], doy, 1)[3]
# }
# 
# dls_2=data.frame(doy=c(1:365), photoperiod=NA)
# for (doy in c(1:365)){
#   dls_2$photoperiod[doy]=insol::daylength(lat=E$lon[1], lon=7.066667, doy, 1)[3]
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
calc_chill=function(phen_data, clim_data, return_all_doy=FALSE){
  #df$X[duplicated(df$X)]
  df=dplyr::full_join(phen_data, clim_data, by=c('site'='site'))
  if (!return_all_doy){
    df=df%>%
    filter((doy<=doy_event&year==year_of_event)|
             (doy>304&year==(year_of_event-1)))
  }else{
    df=df%>%filter((year==year_of_event&doy<=304)|
      doy>304&year==(year_of_event-1))
  }
  df=df%>%
    mutate(
      chill=ifelse(-5<ERA_av_air_temp_C&ERA_av_air_temp_C<5, 1, 0))%>%
    ungroup()%>%
    arrange(Species, site, Biome, year_of_event, tree, date)%>%
    group_by(X, Species, site, Biome, year_of_event, tree)%>%
    mutate(Chill_5_5=cumsum(chill))%>%
    ungroup()%>%
    select(-chill)
  if (!return_all_doy){
    df=df%>%filter(doy_event==doy&year==year_of_event)
  }else{
    df=df%>%
      mutate(doy=ifelse(doy>304, doy-365, doy))
  }
  return(df)
}
#tst=calc_chill(phen_data, clim_data, return_all_doy=TRUE)
#sanity test that you didn't break something with nas or blanks
#sum (duplicated(tst$X))

#define function to accumulate ppt and forcing from Jan 1 to day-of-event
calc_ppt_force=function(phen_data, clim_data, return_all_doy=FALSE){
  df=dplyr::full_join(phen_data, clim_data, by=c('site'='site'))
  if (!return_all_doy){
    df=df%>%filter((doy<=doy_event&!is.na(doy_event)&year_of_event==year))
  }else{
    df=df%>%filter(year_of_event==year)
  }
  df=df%>%
    mutate(
      force=ifelse(ERA_av_air_temp_C>5, 28.4/(1+exp(-.185*(ERA_av_air_temp_C-18.4))), 0))%>%
    ungroup()%>%
    arrange(Species, site, Biome, year_of_event, tree, date)%>%
    group_by(X, Species, site, Biome, year_of_event, tree)%>%
    mutate(Force5s=cumsum(force),
           PrecipDOY=cumsum(total_precipitation))%>%
    ungroup()
  if (!return_all_doy){
  df=df%>%
    filter(doy_event==doy&year==year_of_event)
  }
  return(df)
}
#tst=calc_ppt_force(phen_data, clim_data, return_all_doy=TRUE)


#define function to recalc all the input variables given the climate data
# and some potentially randomized phenology data
est_clim=function(phen_data, clim_data, return_all_doy=FALSE){
  if (!return_all_doy){
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
  }else{
    df=calc_chill(phen_data, clim_data, return_all_doy=TRUE)%>%
      select(-year)%>%
      #calculate forcing and total ppt
      full_join(.,  
                calc_ppt_force(phen_data, clim_data, return_all_doy=TRUE)%>%
                  select(-year)%>%
                  #this will skip photoperiods for the prior yrs Nov data, I think ok?
                  rowwise()%>%
                  mutate(photoperiod=insol::daylength(lat, lon, doy, 1)[3])%>%
                  mutate(photoperiod=ifelse(is.finite(photoperiod), photoperiod, 24)))
  }
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

#generate our best approximation of their climate data using gridded data
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

# -- UNIVARIATE MODEL FITS AND CROSS-VALIDATION  ----

cv_forcing_chilling_photoperiod_site_CV=function(phen_data, clim_data, n_groups=10, prop_drop=.10,
                                                 #need to choose site or year
                                                 by=NULL){
  real_data=est_clim(phen_data, clim_data)
  out=list()
  if (by=='site'){
    sites=unique (real_data$site)
  }
  if (by=='year'){
    years=unique (real_data$year_of_event)
  }
  for (k in 1:n_groups){
    if (by=='site'){
      sites2use=sample(sites, size=round((1-prop_drop)*length(sites)), replace=FALSE)
      fit_data=real_data%>%
        filter(site%in%sites2use)
      val_data=real_data%>%
        filter(!site%in%sites2use)
    }
    if (by=='year'){
      years2use=sample(years, size=round((1-prop_drop)*length(years)), replace=FALSE)
      fit_data=real_data%>%
        filter(year_of_event%in%years2use)
      val_data=real_data%>%
        filter(!year_of_event%in%years2use)
    }
    #same thing but drop one CV
    f1=lme4::lmer(Force5s ~ 1 + (1|Species)+(1|site), data = fit_data,
                  control=lmerControl(optimizer="Nelder_Mead"))
    f2=lme4::lmer(Chill_5_5 ~ 1 + (1|Species)+(1|site), data = fit_data,
                  control=lmerControl(optimizer="Nelder_Mead"))
    f3=lme4::lmer(photoperiod ~ 1 + (1|Species)+(1|site), data = fit_data,
                  control=lmerControl(optimizer="Nelder_Mead"))
    f4=lme4::lmer(doy_event ~ MAT + (1|Species)+(1|site), data = fit_data,
                  control=lmerControl(optimizer="Nelder_Mead"))
    f5=lme4::lmer(doy_event ~ lat + (1|Species)+(1|site), data = fit_data,
                  control=lmerControl(optimizer="Nelder_Mead"))
    f6=lme4::lmer(Force5s ~ lat + (1|Species)+(1|site), data = fit_data,
                  control=lmerControl(optimizer="Nelder_Mead"))
    f7=lme4::lmer(doy_event ~ 1 + (1|Species)+(1|site), data = fit_data,
                  control=lmerControl(optimizer="Nelder_Mead"))
    val_data$preds_force=predict(f1, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_chill=predict(f2, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_photoperiod=predict(f3, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_MAT=predict(f4, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_lat=predict(f5, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_force_lat=predict(f6, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_null=predict(f7, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data=val_data%>%select(X, Species,Type,site,              
                               lat,lon,alt,Biome,              
                               year_of_event,tree,doy_event, preds_force, preds_chill, preds_photoperiod,
                               preds_MAT, preds_lat, preds_force_lat, preds_null)
    # #est Force5s for all doys
    clim_all_doy=est_clim(val_data, clim_data, return_all_doy = TRUE)
    pred_vs_actual_force=clim_all_doy%>%
      rowwise()%>%
      filter(Force5s>preds_force&!is.na(Force5s))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), .groups='drop')%>%
      full_join(clim_all_doy%>%
                  group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
                  summarise(doy_event=unique(doy_event), .groups='drop'))%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 365, doy_pred))
    rmse_force=Metrics::rmse(pred_vs_actual_force$doy_event, pred_vs_actual_force$doy_pred)
    pred_vs_actual_chill=clim_all_doy%>%
      rowwise()%>%
      filter(Chill_5_5>preds_chill&!is.na(Chill_5_5))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), .groups='drop')%>%
      full_join(clim_all_doy%>%
                  group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
                  summarise(doy_event=unique(doy_event), .groups='drop'))%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 303, doy_pred))
    rmse_chill=Metrics::rmse(pred_vs_actual_chill$doy_event, pred_vs_actual_chill$doy_pred)
    pred_vs_actual_photoperiod=clim_all_doy%>%
      rowwise()%>%
      filter(photoperiod>preds_photoperiod&doy>0&doy<=173&!is.na(photoperiod))%>%# take out the late fall vals
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), .groups='drop')%>%
      full_join(clim_all_doy%>%
                   group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
                  summarise(doy_event=unique(doy_event), .groups='drop'))%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 173, doy_pred))
    pred_vs_actual_photoperiod_2=clim_all_doy%>%
      rowwise()%>%
      filter(photoperiod<preds_photoperiod&doy>173&!is.na(photoperiod))%>%# take out the late fall vals
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred_2=min(doy), .groups='drop')%>%
      full_join(clim_all_doy%>%
                  group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
                  summarise(doy_event=unique(doy_event), .groups='drop'))%>%
      mutate(doy_pred_2=ifelse(is.na(doy_pred_2), 173, doy_pred_2))
    pred_vs_actual_photoperiod=full_join(pred_vs_actual_photoperiod, pred_vs_actual_photoperiod_2)%>%
      #pick the nearest one spring or fall
      mutate(doy_pred_all=ifelse(abs(doy_event-doy_pred)<=abs(doy_event-doy_pred_2), doy_pred, doy_pred_2))
    rmse_photoperiod_all=Metrics::rmse(pred_vs_actual_photoperiod$doy_event, pred_vs_actual_photoperiod$doy_pred_all)
    rmse_photoperiod_spring=Metrics::rmse(pred_vs_actual_photoperiod$doy_event, pred_vs_actual_photoperiod$doy_pred)
    rmse_photoperiod_fall=Metrics::rmse(pred_vs_actual_photoperiod$doy_event, pred_vs_actual_photoperiod$doy_pred_2)
    pred_vs_actual_MAT=clim_all_doy%>%
      select(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree, doy_event, preds_MAT)%>%
      distinct()
    rmse_MAT=Metrics::rmse(pred_vs_actual_MAT$doy_event, pred_vs_actual_MAT$preds_MAT)
    pred_vs_actual_lat=clim_all_doy%>%
      select(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree, doy_event, preds_lat)%>%
      distinct()
    rmse_lat=Metrics::rmse(pred_vs_actual_lat$doy_event, pred_vs_actual_lat$preds_lat)
    pred_vs_actual_force_lat=clim_all_doy%>%
      rowwise()%>%
      filter(Force5s>preds_force_lat&!is.na(Force5s))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), doy_event=unique(doy_event), .groups='drop')%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 365, doy_pred))
    rmse_force_lat=Metrics::rmse(pred_vs_actual_force_lat$doy_event, pred_vs_actual_force_lat$doy_pred)
    pred_vs_actual_null=clim_all_doy%>%
      select(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree, doy_event, preds_null)%>%
      distinct()
    rmse_null=Metrics::rmse(pred_vs_actual_null$doy_event, pred_vs_actual_null$preds_null)
    
    out[[k]]=list(sites2use=ifelse(by=='site', paste0(sites2use, collapse='|'), NA), 
                  years2use=ifelse(by=='year', paste0(years2use, collapse='|'), NA), rmse_force=rmse_force, rmse_chill=rmse_chill,
                  rmse_photoperiod_all=rmse_photoperiod_all, rmse_photoperiod_spring=rmse_photoperiod_spring,
                  rmse_photoperiod_fall=rmse_photoperiod_fall,
                  rmse_MAT=rmse_MAT, rmse_lat=rmse_lat,
                  rmse_force_lat=rmse_force_lat,
                  rmse_null=rmse_null)
  }
  return(out)
}

#run - takes a while
tst_site=cv_forcing_chilling_photoperiod_site_CV(phen_data, clim_data, n_groups=100, prop_drop=0.1,
                                            by='site')

tst_year=cv_forcing_chilling_photoperiod_site_CV(phen_data, clim_data, n_groups=100, prop_drop=0.1,
                                            by='year')

#convert to df to plot
#there's probably better list indexing way to do this
# it just make a ton of duplicates that you can then ditch

#this for the boxplots
tst_site_df=data.table::rbindlist(tst_site, idcol='group')%>%
  select(-sites2use, -years2use)%>%
  distinct()%>%
  pivot_longer(cols=-group)%>%
  rename(model=name, rmse=value)%>%
  mutate(model=gsub('rmse_', '', model))

#this for the density plots
tst_site_df_wide=data.table::rbindlist(tst_site, idcol='group')%>%
  select(-years2use)%>%
  rename(null_val=rmse_null)%>%
  distinct()%>%  mutate_at(vars(contains('rmse')), .funs = rlang::list2(~. - null_val))%>%
  select(-null_val)%>%
  pivot_longer(cols=c(-group, -sites2use))%>%
  rename(model=name, rmse_delta=value)%>%
  mutate(model=gsub('rmse_', '', model))

tst_year_df=data.table::rbindlist(tst_year, idcol='group')%>%
  select(-years2use, -sites2use)%>%
  distinct()%>%
  pivot_longer(cols=-group)%>%
  rename(model=name, rmse=value)%>%
  mutate(model=gsub('rmse_', '', model))

tst_year_df_wide=data.table::rbindlist(tst_year, idcol='group')%>%
  select(-sites2use)%>%
  rename(null_val=rmse_null)%>%
  distinct()%>%
  mutate_at(vars(contains('rmse')), .funs = rlang::list2(~. - null_val))%>%
  select(-null_val)%>%
  pivot_longer(cols=c(-group, -years2use))%>%
  rename(model=name, rmse_delta=value)%>%
  mutate(model=gsub('rmse_', '', model))

#note the chill and fall photoperiod so bad you can't even get there
null_year_model_plot = ggplot(data=tst_year_df_wide)+
                      geom_density(aes(x=rmse_delta, group=model, fill=model), 
                        alpha=0.5, adjust=2)+
  geom_vline(xintercept=0, linetype='dotted')+
  #facet_grid(~model)+xlim(-10,10)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(null_year_model_plot,
       file='plots/null_year_model_plot.jpg',
       width=10, height=6)

null_year_model_plot_full = ggplot(data=tst_year_df_wide)+
  geom_density(aes(x=rmse_delta, group=model, fill=model), 
               alpha=0.5, adjust=2)+
  geom_vline(xintercept=0, linetype='dotted')+
  facet_grid(~model)+#xlim(-10,10)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(null_year_model_plot_full,
       file='plots/null_year_model_plot_full.jpg',
       width=10, height=6)


null_site_model_plot = ggplot(data=tst_site_df_wide)+
  geom_density(aes(x=rmse_delta, group=model, fill=model), 
               alpha=0.5, adjust=2)+
  geom_vline(xintercept=0, linetype='dotted')+
  facet_grid(~model)+xlim(-10,10)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(null_site_model_plot,
       file='plots/null_site_model_plot.jpg',
       width=10, height=6)

#or vis like this?
ggplot(data=tst_year_df, aes(color=model, x=model, y=rmse))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data=tst_site_df, aes(color=model, x=model, y=rmse))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

