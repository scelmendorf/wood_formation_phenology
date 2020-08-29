## 28 Aug 2020 ##
## SCE ##

# -- NOTES ON DATA ----
#1
# Data and metadata do not make clear what the units of their ppt are.
# Mean annual is ~908 so guessing its mm since "globally averaged annual
# precipitation over land is 715 mm (by wikipedia). Ppt totals are calculated
# here with that unit conversion assumption but not used in the final 
# analyses since it was not a major thrust of the paper


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
## working directory should be root of this project, all paths then relative
rm(list=ls())
lapply(paste("package:", names(sessionInfo()$otherPkgs), sep=""), detach, character.only = TRUE, unload = TRUE)

set.seed (123456)

#path to the reanalysis met data
#if you do not wish to reextract using GEE  (takes a few hours and requires an
#account, you can download from gdrive here:
#https://drive.google.com/file/d/165Owv5gtYrBYTmPJXwPGJx9TSAGP8iLo/view?usp=sharing)

#this file is too big to be archived on git so you will need to manually
# add it to this location for the code to run
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

#this is used as the input to google earth engine code
#see code scripts/GEE_ERA_extraction.txt
write.csv(E_2_earth_engine, 'data/E_2_earth_engine.csv',
          row.names=FALSE)
}

# -- READ IN AND FORMAT THE CLIM DATA  ----

ERA_clim=read.csv(ERA_clim_path)

#calculate MAT off the ERA data
# use that as a bias correction for ERAs daily temp data
# following the logic of 
# (1) Bias correct using the same year if available
# (2) Bias correct using the next year if same year not available
# (3) Bias correct using the average of all years if same/next not available
# always by site location. For temp, adjusted based on the difference in means
# for ppt, adjusted based on the ratio of differences as ppt cannot be negative

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

#note there are some oddities in S1 data file, e.g.
# looks like the ppt for all some trees trees in this site*yr
# is missing but since it is a site-year level variable and constant
# for the rest of the measurements at that site*year, assume it was just
# inadvertently deleted? Does not really affect the results
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
# or all days of year (if return_all_doy=TRUE)
calc_chill=function(phen_data, clim_data, return_all_doy=FALSE){
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

#sanity test
#tst=calc_chill(phen_data, clim_data, return_all_doy=TRUE)
#sum (duplicated(tst$X))

#define function to accumulate ppt and forcing from Jan 1 to day-of-event
# or all days of year (if return_all_doy=TRUE)
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
# and phenology data
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
    mutate(month_prior_event=ifelse(day_of_month<=15,
                                    lubridate::month(date_event)-1, lubridate::month(date_event)))%>%
    rowwise()%>%
    mutate(jday_event=insol::JDymd(lubridate::year(date_event), lubridate::month(date_event), day_of_month))%>%
    mutate(photoperiod=insol::daylength(lat, lon, jday_event, 1)[3])%>%
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
    # surely there is a better way to deal with the type conversions here
    mutate(scPDSI_p=ifelse(scPDSI_p==-9999, NA, scPDSI_p))#%>%
  }else{
    df=calc_chill(phen_data, clim_data, return_all_doy=TRUE)%>%
      select(-year)%>%
      #calculate forcing and total ppt
      full_join(.,  
                calc_ppt_force(phen_data, clim_data, return_all_doy=TRUE)%>%
                  select(-year)%>%
                  #this will skip photoperiods for the prior yrs Nov data
                  # ok since that is never even close to phenodates
                  rowwise()%>%
                  mutate(photoperiod=insol::daylength(lat, lon, doy, 1)[3])%>%
                  mutate(photoperiod=ifelse(is.finite(photoperiod), photoperiod, 24)))
  }
  return(df)
}

#make single climate file with all adjusted clim variables
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

phen_data=E%>%rename(year_of_event=year, doy_event=DOY)%>%
  #delete the climate variables used in the manuscript; we will
  # replace with those derived from ERA for consistency
select(-Chill0_5, -Chill_5_5,    
       -Chill_5_0, -Chill_10_0, -Force5s, -photoperiod, -PrecipDOY, -scPDSI_p)


clim_data=ERA_clim%>%select(site, year, date, ERA_av_air_temp_C, total_precipitation)%>%
  mutate(doy=lubridate::yday(date))

#generate our best approximation of their climate data using gridded data
# i.e. this is real phenology data, but with the forcing variables calculated
# from ERA rather than local met data
real_data=est_clim(phen_data, clim_data)

#sanity checks that we calculated right
#our guesstimates of chilling/forcingpretty close not exact but not super biased
valid_ERA_plot_force=ggplot(data=data.frame(observed_GDD=real_data$Force5s[order(real_data$X)],
           estimated_GDD=E$Force5s[order(E$X)]),
       aes(x=estimated_GDD, y=observed_GDD))+
         geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme(panel.background = element_blank())+
  labs(x = 'estimated GDD from gridded data', y='observed GDD from local sources')
  
ggsave("plots/valid_ERA_plot_force.jpg", valid_ERA_plot_force, width=5, height=5)

#rmse for reference
#Metrics::rmse(real_data$Force5s[order(real_data$X)], E$Force5s[order(E$X)])

valid_ERA_plot_chill=ggplot(data=data.frame(observed_Chill=real_data$Chill_5_5[order(real_data$X)],
                                      estimated_Chill=E$Chill_5_5[order(E$X)]),
                      aes(x=estimated_Chill, y=observed_Chill))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme(panel.background = element_blank())+
  labs(x = 'estimated Chill from gridded data', y='observed Chill from local sources')
ggsave("plots/valid_ERA_plot_chill.jpg", valid_ERA_plot_chill, width=5, height=5)

#check we have photoperiod calculated correctly
valid_photoperiod_plot=ggplot(data=data.frame(S1_photoperiod=real_data$photoperiod[order(real_data$X)],
                                            estimated_photoperiod=E$photoperiod[order(E$X)]),
                            aes(x=estimated_photoperiod, y=S1_photoperiod))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme(panel.background = element_blank())+
  labs(x = 'estimated photoperiod from INSOL package', y='photoperiod from table S1')
ggsave("plots/valid_ERA_plot_photoperiod.jpg", valid_ERA_plot_chill, width=5, height=5)

#these should be exact match bc I just looked for the right month
valid_scPDSI_p_plot=ggplot(data=data.frame(observed_scPDSI_p=real_data$scPDSI_p[order(real_data$X)],
                                              estimated_scPDSI_p=E$scPDSI_p[order(E$X)]),
                              aes(x=estimated_scPDSI_p, y=observed_scPDSI_p))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme(panel.background = element_blank())+
  labs(x = 'estimated _scPDSI_p', y='observed_scPDSI_p')

# -- UNIVARIATE MODEL FITS AND CROSS-VALIDATION  ----

cv_forcing_chilling_photoperiod_site_CV=function(phen_data, clim_data, n_groups=10, prop_drop=.10,
                                                 #need to choose site or year
                                                 by=NULL, choose=FALSE){
  #' This function does iterative blocked cross validation of phenology predictions
  #' @param phen_data phenology data
  #' @param clim_data daily climate data for all sites in phen_data
  #' @param n_groups number of cv iteractions to do (set only if not using choose method)
  #' @param prop_drop percentage of blocks to drop in each iteration
  #' @param by either site or year - how to block subsampling for cross validation
  #' @param choose if TRUE, do all combinations. Caution this could run for a long time
  real_data=est_clim(phen_data, clim_data)
  out=list()
  if (by=='site'){
    sites=unique (real_data$site)
    #don't do all combinations of site, it will take forever
  }
  if (by=='year'){
    years=unique (real_data$year_of_event)
  }
  if (choose){
    all_combinations=combn(c(1:length(years)),round((1-prop_drop)*length(years)))
    n_groups=ncol(all_combinations)
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
      if (choose){
        years2use=years[all_combinations[,k]]
      }else{
      years2use=sample(years, size=round((1-prop_drop)*length(years)), replace=FALSE)
      }
      fit_data=real_data%>%
        filter(year_of_event%in%years2use)
      val_data=real_data%>%
        filter(!year_of_event%in%years2use)
    }
    #same thing but drop one CV
    f1=lme4::lmer(Force5s ~ 1 + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f2=lme4::lmer(Chill_5_5 ~ 1 + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f3=lme4::lmer(photoperiod ~ 1 + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f4=lme4::lmer(doy_event ~ MAT + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f5=lme4::lmer(doy_event ~ lat + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f6=lme4::lmer(Force5s ~ lat + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f7=lme4::lmer(doy_event ~ 1 + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
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

#run cross validations - takes a while
cv_site=cv_forcing_chilling_photoperiod_site_CV(phen_data, clim_data, n_groups=100, prop_drop=0.1,
                                            by='site')

cv_year=cv_forcing_chilling_photoperiod_site_CV(phen_data, clim_data, prop_drop=0.1,
                                            by='year', choose=TRUE)

#convert to df to plot
#there's surely a faster list indexing way to do this
#but this works
# it just make a ton of duplicates that you can then ditch with distinct()

#this for the density plots
cv_site_df_wide=data.table::rbindlist(cv_site, idcol='group')%>%
  select(-years2use)%>%
  rename(null_val=rmse_null)%>%
  distinct()%>%  mutate_at(vars(contains('rmse')), .funs = rlang::list2(~. - null_val))%>%
  select(-null_val)%>%
  pivot_longer(cols=c(-group, -sites2use))%>%
  rename(model=name, rmse_delta=value)%>%
  mutate(model=gsub('rmse_', '', model))

#report only summary stats for chill as putting them in the figure
# they are so bad the x axis range has to be huge
mean(cv_site_df_wide$rmse_delta[cv_site_df_wide$model=='chill'])
quantile(cv_site_df_wide$rmse_delta[cv_site_df_wide$model=='chill'])

cv_year_df_wide=data.table::rbindlist(cv_year, idcol='group')%>%
  select(-sites2use)%>%
  rename(null_val=rmse_null)%>%
  distinct()%>%
  mutate_at(vars(contains('rmse')), .funs = rlang::list2(~. - null_val))%>%
  select(-null_val)%>%
  pivot_longer(cols=c(-group, -years2use))%>%
  rename(model=name, rmse_delta=value)%>%
  mutate(model=gsub('rmse_', '', model))

#these so bad just report the summaries
mean(cv_year_df_wide$rmse_delta[cv_year_df_wide$model=='chill'])
quantile(cv_year_df_wide$rmse_delta[cv_year_df_wide$model=='chill'])

#can look at percentage that are better or worse than null if desired
cv_year_df_wide%>%
  group_by(model)%>%
  summarise(sum(rmse_delta>0)/dplyr::n())

cv_site_df_wide%>%
  group_by(model)%>%
  summarise(sum(rmse_delta>0)/dplyr::n())

#note the chill and fall photoperiod so bad you can't even get there
null_year_model_plot = ggplot(data=cv_year_df_wide%>%
                                filter(!model%in%c('chill', 'photoperiod_fall',
                                                   'force_lat', 'photoperiod_spring'))%>%
                                mutate(
                                  model=case_when(
                                    model=='force' ~ 'Forcing',
                                    model=='lat' ~ 'Latitude',
                                    model=='photoperiod_all' ~ 'Photoperiod',
                                    model=='MAT' ~ 'MAT',
                                    TRUE ~ 'NA')
                                ))+
                      geom_density(aes(x=rmse_delta, group=model, fill=model, color=model), 
                        alpha=0.5, adjust=2)+
  theme(panel.background = element_blank())+
  geom_vline(xintercept=0, linetype='dotted')+
  facet_grid(rows = 'model')+
  #xlim(-15,15)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = expression(delta[rmse(model - null)]))+
  guides(fill=FALSE, color=FALSE)

#add intuitive label axis
temporal_cv=gridExtra::grid.arrange(
  null_year_model_plot,
  nrow = 1,
  top = "Predictive performace (interannual variation)",
  bottom = textGrob(
    expression(better %<->% worse),
    gp = gpar(fontface = 3, fontsize = 25),
    just='center'
  ))

ggsave("plots/temporal_cv.jpg", temporal_cv)

#repeat for site
null_site_model_plot = ggplot(data=cv_site_df_wide%>%
                                filter(!model%in%c('chill', 'photoperiod_fall',
                                                   'force_lat', 'photoperiod_spring'))%>%
                                mutate(
                                  model=case_when(
                                    model=='force' ~ 'Forcing',
                                    model=='lat' ~ 'Latitude',
                                    model=='photoperiod_all' ~ 'Photoperiod',
                                    model=='MAT' ~ 'MAT',
                                    TRUE ~ 'NA')
                                ))+
  geom_density(aes(x=rmse_delta, group=model, fill=model, color=model), 
               alpha=0.5, adjust=2)+
  theme(panel.background = element_blank())+
  geom_vline(xintercept=0, linetype='dotted')+
  facet_grid(rows = 'model')+
  #xlim(-15,15)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = expression(delta[rmse(model - null)]))+
  guides(fill=FALSE, color=FALSE)

#add intuitive label axis
spatial_cv=gridExtra::grid.arrange(
  null_site_model_plot,
  nrow = 1,
  top = "Predictive performace (spatial variation)",
  bottom = textGrob(
    expression(better %<->% worse),
    gp = gpar(fontface = 3, fontsize = 25),
    #hjust = 1,
    just='center'
    #x = 0.1
  ))

ggsave("plots/spatial_cv.jpg", spatial_cv)

#figure out range of year models so you can glue them together
#into a single multi-panel figure but reducing
# labeling redundancy in cowplot
ylims=ggplot_build(null_year_model_plot)$layout$panel_scales_y[[1]]$range$range

spatial_cv_2panel=null_site_model_plot+
 theme(
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )+ylim(ylims)+
  labs(title='spatial variation')

spatial_cv_2panel=gridExtra::grid.arrange(
  spatial_cv_2panel,
  nrow = 1,
  bottom = textGrob(
    expression(better %<->% worse),
    gp = gpar(fontface = 3, fontsize = 25),
    just='center'
  ))

temporal_cv_2panel=null_year_model_plot+
  labs(y=NULL)+
  labs(title='temporal variation')+
  theme(axis.ticks = element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5))

temporal_cv_2panel=gridExtra::grid.arrange(
  temporal_cv_2panel,
  nrow = 1,
  #top = "Predictive performace (spatial variation)",
  bottom = textGrob(
    expression(better %<->% worse),
    gp = gpar(fontface = 3, fontsize = 25),
    just='center'
  ))

panelfig_2=cowplot::plot_grid(spatial_cv_2panel, temporal_cv_2panel)

title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "Model predictive performance",
    fontface = 'bold',
    hjust = 0.5
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

panelfig_2=cowplot::plot_grid(
  title, panelfig_2,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

ggsave("plots/spatial_temporal_cv.jpg", panelfig_2, width=6, height=5)

