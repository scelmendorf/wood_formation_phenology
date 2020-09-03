## 2 Sept 2020 ##
## SCE ##

#sce fix the chill_5_5 april 30 has a ton of NAs in ther
# maybe because it's often 0 not usre

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

#note dropping PPT for now as not a major part of the story
#site, year, date, GDD to date, chill to date, and chill to april 30

#sanity test for still have 365 d in this year
#View (ERA_clim%>%filter(site=='5T1'&year==2001))

clim_data=ERA_clim%>%select(site, year, date, ERA_av_air_temp_C, total_precipitation)%>%
  mutate(doy=lubridate::yday(date))%>%
  #add in latlong
  left_join(., E%>%select(site, lat, lon)%>%distinct())

#cleanup
rm(ERA_clim)

#define function to accumulate chill from Nov to day-of-event
#force from Jan 1 to day of event
# chill from Nov thru April 30
#photoperiod for each doy


calc_chill_force_p_all=function(clim_data){
  df_chill=clim_data%>%
    mutate(doy=lubridate::yday(date))%>%
    mutate(water_year=ifelse(doy>304, year+1, year))%>%
    mutate(
      chill=ifelse(-5<ERA_av_air_temp_C&ERA_av_air_temp_C<5, 1, 0))%>%
    ungroup()%>%
    arrange(site, date)%>%
    group_by(site, water_year)%>%
    mutate(Chill_5_5=cumsum(chill), nday=dplyr::n())%>%
    #only do where there are at least 365 d in the year
    mutate(Chill_5_5=ifelse(nday<365, NA, Chill_5_5))%>%
    #do not entertain chill sums that occur in Nov-dec
    mutate(Chill_5_5=ifelse(doy>304, NA, Chill_5_5))%>%
    select(water_year, site, Chill_5_5, doy)%>%
    rename(year=water_year)
  #add in the GDD sums
  df_force=clim_data%>%
    mutate(
      force=ifelse(ERA_av_air_temp_C>5, 28.4/(1+exp(-.185*(ERA_av_air_temp_C-18.4))), 0))%>%
    ungroup()%>%
    arrange(site, date)%>%
    group_by(site, year)%>%
    mutate(Force5s=cumsum(force),
           nday=dplyr::n())%>%
    ungroup()%>%
    mutate(doy=lubridate::yday(date))%>%
    select(site, year, date, Force5s, doy)
  #add in photoperiod
  df_photoperiod=clim_data%>%
    rowwise()%>%
    mutate(jday=insol::JDymd(lubridate::year(date), lubridate::month(date), 
                                    lubridate::day(date)))%>%
    rowwise()%>%
    #warnings on 24 hour photoperiod ok, solved in next line
    mutate(photoperiod=suppressWarnings(insol::daylength(lat, lon, jday, 1)[3]))%>%
    mutate(photoperiod=ifelse(is.finite(photoperiod), photoperiod, 24))%>%
    ungroup()%>%
    mutate(doy=lubridate::yday(date))%>%
    select(site, date, doy, photoperiod)
  
  #combine  
  df_all=df_force%>%
    left_join(., df_chill)%>%
    left_join(., df_photoperiod)%>%
    #add in the April 30 chill
    left_join(., 
              df_chill%>%
              filter(doy==121)%>%
              rename(Chill_5_5_april_30=Chill_5_5)%>%
              select(-doy)
              )
  return(df_all)
}



#run in single DOY model for the clim_data we have
clim_data=calc_chill_force_p_all(
  clim_data
)

#define function to either return clim up to a DOY event, or over all doy
#as needed
est_clim=function(phen_data, clim_data, return_all_doy=FALSE){
  if (return_all_doy){
    df=full_join(phen_data, clim_data, by=c('year_of_event'='year', 'site'='site'))
  }else{
    df=left_join(phen_data, 
                 clim_data%>%rename(doy_event=doy),
                 by=c('year_of_event'='year', 'site'='site', 'doy_event'='doy_event'))
  }
}

# -- READ IN THE REAL DATA AND MAKE CLIMATE DATA THAT MATCHES IT  ----

phen_data=E%>%rename(year_of_event=year, doy_event=DOY)%>%
  #delete the climate variables used in the manuscript; we will
  # replace with those derived from ERA for consistency
  select(-Chill0_5, -Chill_5_5,    
         -Chill_5_0, -Chill_10_0, -Force5s, -photoperiod, -PrecipDOY, -scPDSI_p)

#generate our best approximation of their climate data using gridded data
# i.e. this is real phenology data, but with the forcing variables calculated
# from ERA rather than local met data
real_data=est_clim(phen_data, clim_data, return_all_doy=FALSE)


#sanity checks that using ERA clim is more or less accurate
# and our photoperiod calcs are correct 
#our guesstimates of chilling/forcingpretty close not exact but not super biased
valid_ERA_plot_force=ggplot(data=data.frame(observed_Force5s=real_data$Force5s[order(real_data$X)],
           estimated_Force5s=E$Force5s[order(E$X)]),
       aes(x=estimated_Force5s, y=observed_Force5s))+
         geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme(panel.background = element_rect(fill = "transparent"))+
  labs(x = 'estimated forcing using gridded data', y='observed forcing from local sources')
  
ggsave("plots/valid_ERA_plot_force.jpg", valid_ERA_plot_force, width=5, height=5)

valid_ERA_plot_chill=ggplot(data=data.frame(observed_Chill_5_5=real_data$Chill_5_5[order(real_data$X)],
                                      estimated_Chill_5_5=E$Chill_5_5[order(E$X)]),
                      aes(x=estimated_Chill_5_5, y=observed_Chill_5_5))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme(panel.background = element_blank())+
  labs(x = 'estimated chilling from gridded data', y='observed chilling  from local sources')
ggsave("plots/valid_ERA_plot_chill.jpg", valid_ERA_plot_chill, width=5, height=5)

#check we have photoperiod calculated correctly
valid_photoperiod_plot=ggplot(data=data.frame(observed_photoperiod=real_data$photoperiod[order(real_data$X)],
                                            estimated_photoperiod=E$photoperiod[order(E$X)]),
                            aes(x=estimated_photoperiod, y=observed_photoperiod))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme(panel.background = element_blank())+
  labs(x = 'estimated photoperiod from INSOL package', y='photoperiod from table S1')
ggsave("plots/valid_ERA_plot_photoperiod.jpg", valid_photoperiod_plot, width=5, height=5)



# -- UNIVARIATE MODEL FITS AND CROSS-VALIDATION  ----
#define cross-validation function
cv_forcing_chilling_photoperiod_CV=function(phen_data, clim_data, n_groups=10, prop_drop=.10,
                                                 #need to choose by of site or year or year_within_site
                                                 # year within site seemed more logical than year so 
                                                 # used that in final sims
                                                 by=NULL, choose=FALSE){
  #' This function does iterative blocked cross validation of phenology predictions
  #' @param phen_data phenology data
  #' @param clim_data daily climate data for all sites in phen_data
  #' @param n_groups number of cv iteractions to do (set only if not using choose method)
  #' @param prop_drop percentage of blocks to drop in each iteration, ignored if choose=TRUE
  #' @param by either site or year or year_within_site - how to block subsampling for cross validation
  #' @param choose if TRUE, do all combinations. Caution this could run for a long time
  real_data=est_clim(phen_data, clim_data, return_all_doy=FALSE)
  #centering the chill sums to april 30 to help with convergence
  real_data$Chill_5_5_april_30_scale=scale(real_data$Chill_5_5_april)
  real_data$lat_scale=scale(real_data$lat)
  real_data$mat_scale=scale(real_data$MAT)
  out=list()
  if (by=='site'){
    sites=unique (real_data$site)
    #don't do all combinations of site, it will take forever
  }
  if (by=='year'){
    years=unique (real_data$year_of_event)
  }
  if (by=='year_within_site'){
    yearsites=real_data%>%select(year_of_event, site)%>%distinct()
  }
  
  if (choose){
    all_combinations=combn(c(1:length(years)),round((1-prop_drop)*length(years)))
    n_groups=ncol(all_combinations)
  }
  
  for (k in 1:n_groups){
    pb <- txtProgressBar(min = 0, max = n_groups, style = 3)
    setTxtProgressBar(pb, k)
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
    
    if (by=='year_within_site'){
      val_data=yearsites%>%
        group_by(site)%>%
        summarise(year_of_event=sample(year_of_event, 1),
                  nyr=dplyr::n())%>%
        filter(nyr>1)%>%
        select(-nyr)
      fit_data=real_data%>%
        anti_join(., val_data)
      val_data=real_data%>%
        inner_join(., val_data)
    }
    
    #fit a variety of models to the train subset
    f1=lme4::lmer(Force5s ~ 1 + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f1a=lme4::lmer(Force5s ~ lat + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f1b=lme4::lmer(Force5s ~ Chill_5_5_april_30_scale + (1|Species)+(1|site), data = fit_data,
                   control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f1c=lme4::lmer(Force5s ~ Chill_5_5_april_30_scale + (Chill_5_5_april_30_scale|Species)+(1|site), data = fit_data,
                   control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f1d=lme4::lmer(Force5s ~ Chill_5_5_april_30_scale + (1|Species)+(Chill_5_5_april_30_scale|site), data = fit_data,
                   control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f1e=lme4::lmer(Force5s ~ Chill_5_5_april_30_scale + lat + (1|Species)+(1|site), data = fit_data,
                   control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f1f=lme4::lmer(Force5s ~ Chill_5_5_april_30_scale * lat + (1|Species)+(1|site), data = fit_data,
                   control=lme4::lmerControl(optimizer="Nelder_Mead"))
    
    f2=lme4::lmer(Chill_5_5 ~ 1 + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f3=lme4::lmer(photoperiod ~ 1 + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f3a=lme4::lmer(photoperiod ~ MAT + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    
    #this is new
    f3a1=lme4::lmer(photoperiod ~ mat_scale + (mat_scale|Species)+(1|site), data = fit_data,
                   control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f3b=lme4::lmer(photoperiod ~ lat + (1|Species)+(1|site), data = fit_data,
                   control=lme4::lmerControl(optimizer="Nelder_Mead"))

    #this one was breaking when not scaled
    f3c=lme4::lmer(photoperiod ~ lat_scale + (lat_scale|Species)+(1|site), data = fit_data,
                   control=lme4::lmerControl(optimizer="Nelder_Mead"))
    f4=lme4::lmer(doy_event ~ MAT + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    #this is new
    f4a=lme4::lmer(doy_event ~ mat_scale + (mat_scale|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    
    f5=lme4::lmer(doy_event ~ lat + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    #new
    f5a= lme4::lmer(doy_event ~ lat_scale + (lat_scale|Species)+(1|site), data = fit_data,
               control=lme4::lmerControl(optimizer="Nelder_Mead"))

    f6=lme4::lmer(doy_event ~ 1 + (1|Species)+(1|site), data = fit_data,
                  control=lme4::lmerControl(optimizer="Nelder_Mead"))
    
    #generate predictions on variety of scales for test subset
    #forcing suite
    val_data$preds_force=predict(f1, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_forcea=predict(f1a, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_forceb=predict(f1b, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_forcec=predict(f1c, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_forced=predict(f1d, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_forcee=predict(f1e, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_forcef=predict(f1f, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    
    #chill
    val_data$preds_chill=predict(f2, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    
    #photoperiod suite
    val_data$preds_photoperiod=predict(f3, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_photoperioda=predict(f3a, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_photoperioda1=predict(f3a1, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    
    val_data$preds_photoperiodb=predict(f3b, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_photoperiodc=predict(f3c, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    
    #model doy directly
    val_data$preds_MAT=predict(f4, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_MATa=predict(f4a, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    
    val_data$preds_lat=predict(f5, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data$preds_lata=predict(f5a, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    
    
    val_data$preds_null=predict(f6, newdata = val_data, allow.new.levels = TRUE, re.form = NULL)
    val_data=val_data%>%select(X, Species,Type,site,              
                               lat,lon,alt,Biome,              
                               year_of_event,tree,doy_event, preds_force, preds_forcea,
                               preds_forceb,preds_forcec,preds_forced,preds_forcee,preds_forcef,
                               preds_chill, preds_photoperiod,
                               preds_photoperioda,preds_photoperioda1,preds_photoperiodb,preds_photoperiodc,
                               preds_MAT, preds_MATa, preds_lat, preds_lata, preds_null)
    
    #daily estimate clim for val data subset
    clim_all_doy=est_clim(val_data, clim_data, return_all_doy = TRUE)
    
    pred_vs_actual_force=clim_all_doy%>%
      rowwise()%>%
      filter(Force5s>preds_force&!is.na(Force5s))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), .groups='drop')%>%
      #rejoin to real data
      full_join(val_data%>%
                  select(!matches('preds')))%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 365, doy_pred))
    
    rmse_force=Metrics::rmse(pred_vs_actual_force$doy_event, pred_vs_actual_force$doy_pred)
    
    pred_vs_actual_chill=clim_all_doy%>%
      rowwise()%>%
      filter(Chill_5_5>preds_chill&!is.na(Chill_5_5))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), .groups='drop')%>%
      full_join(val_data%>%
                  select(!matches('preds')))%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 303, doy_pred))
    
    rmse_chill=Metrics::rmse(pred_vs_actual_chill$doy_event, pred_vs_actual_chill$doy_pred)
    
    pred_vs_actual_photoperiod=clim_all_doy%>%
      rowwise()%>%
      filter(photoperiod>preds_photoperiod&doy>0&doy<=173&!is.na(photoperiod))%>%# take out the late fall vals
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), .groups='drop')%>%
      full_join(val_data%>%
                  select(!matches('preds')))%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 173, doy_pred))
    
    pred_vs_actual_photoperiod_2=clim_all_doy%>%
      rowwise()%>%
      filter(photoperiod<preds_photoperiod&doy>173&!is.na(photoperiod))%>%# include only fall vals
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred_2=min(doy), .groups='drop')%>%
      full_join(val_data%>%
                  select(!matches('preds')))%>%
      mutate(doy_pred_2=ifelse(is.na(doy_pred_2), 173, doy_pred_2))
    pred_vs_actual_photoperiod=full_join(pred_vs_actual_photoperiod, pred_vs_actual_photoperiod_2)%>%
      #pick the nearest one spring or fall
      mutate(doy_pred_all=ifelse(abs(doy_event-doy_pred)<=abs(doy_event-doy_pred_2), doy_pred, doy_pred_2))
    rmse_photoperiod_all=Metrics::rmse(pred_vs_actual_photoperiod$doy_event, pred_vs_actual_photoperiod$doy_pred_all)
    rmse_photoperiod_spring=Metrics::rmse(pred_vs_actual_photoperiod$doy_event, pred_vs_actual_photoperiod$doy_pred)
    rmse_photoperiod_fall=Metrics::rmse(pred_vs_actual_photoperiod$doy_event, pred_vs_actual_photoperiod$doy_pred_2)
    
    #lat and photoperiod
    pred_vs_actual_photoperioda=clim_all_doy%>%
      rowwise()%>%
      filter(photoperiod>preds_photoperioda&doy>0&doy<=173&!is.na(photoperiod))%>%#
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), .groups='drop')%>%
      full_join(val_data%>%
                  select(!matches('preds')))%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 173, doy_pred))
    pred_vs_actual_photoperiod_2a=clim_all_doy%>%
      rowwise()%>%
      filter(photoperiod<preds_photoperioda&doy>173&!is.na(photoperiod))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred_2=min(doy), .groups='drop')%>%
      full_join(val_data%>%
                  select(!matches('preds')))%>%
      mutate(doy_pred_2=ifelse(is.na(doy_pred_2), 173, doy_pred_2))
    pred_vs_actual_photoperioda=full_join(pred_vs_actual_photoperioda, pred_vs_actual_photoperiod_2a)%>%
      #pick the nearest one spring or fall
      mutate(doy_pred_all=ifelse(abs(doy_event-doy_pred)<=abs(doy_event-doy_pred_2), doy_pred, doy_pred_2))
    rmse_photoperiod_alla=Metrics::rmse(pred_vs_actual_photoperioda$doy_event, pred_vs_actual_photoperioda$doy_pred_all)
    rmse_photoperiod_springa=Metrics::rmse(pred_vs_actual_photoperioda$doy_event, pred_vs_actual_photoperioda$doy_pred)
    rmse_photoperiod_falla=Metrics::rmse(pred_vs_actual_photoperioda$doy_event, pred_vs_actual_photoperioda$doy_pred_2)
    
    
    ##
    #photoperiod a1
    pred_vs_actual_photoperioda1=clim_all_doy%>%
      rowwise()%>%
      filter(photoperiod>preds_photoperioda1&doy>0&doy<=173&!is.na(photoperiod))%>%#
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), .groups='drop')%>%
      full_join(val_data%>%
                  select(!matches('preds')))%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 173, doy_pred))
    pred_vs_actual_photoperiod_2a1=clim_all_doy%>%
      rowwise()%>%
      filter(photoperiod<preds_photoperioda1&doy>173&!is.na(photoperiod))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred_2=min(doy), .groups='drop')%>%
      full_join(val_data%>%
                  select(!matches('preds')))%>%
      mutate(doy_pred_2=ifelse(is.na(doy_pred_2), 173, doy_pred_2))
    pred_vs_actual_photoperioda1=full_join(pred_vs_actual_photoperioda1, pred_vs_actual_photoperiod_2a1)%>%
      #pick the nearest one spring or fall
      mutate(doy_pred_all=ifelse(abs(doy_event-doy_pred)<=abs(doy_event-doy_pred_2), doy_pred, doy_pred_2))
    rmse_photoperiod_alla1=Metrics::rmse(pred_vs_actual_photoperioda1$doy_event, pred_vs_actual_photoperioda1$doy_pred_all)
    rmse_photoperiod_springa1=Metrics::rmse(pred_vs_actual_photoperioda1$doy_event, pred_vs_actual_photoperioda1$doy_pred)
    rmse_photoperiod_falla1=Metrics::rmse(pred_vs_actual_photoperioda1$doy_event, pred_vs_actual_photoperioda1$doy_pred_2)
    
    ##
    
    ##next lat by photoperiod
    pred_vs_actual_photoperiodb=clim_all_doy%>%
      rowwise()%>%
      filter(photoperiod>preds_photoperiodb&doy>0&doy<=173&!is.na(photoperiod))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), .groups='drop')%>%
      full_join(val_data%>%
                  select(!matches('preds')))%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 173, doy_pred))
    pred_vs_actual_photoperiod_2b=clim_all_doy%>%
      rowwise()%>%
      filter(photoperiod<preds_photoperiodb&doy>173&!is.na(photoperiod))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred_2=min(doy), .groups='drop')%>%
      full_join(val_data%>%
                  select(!matches('preds')))%>%
      mutate(doy_pred_2=ifelse(is.na(doy_pred_2), 173, doy_pred_2))
    pred_vs_actual_photoperiodb=full_join(pred_vs_actual_photoperiodb, pred_vs_actual_photoperiod_2b)%>%
      #pick the nearest one spring or fall
      mutate(doy_pred_all=ifelse(abs(doy_event-doy_pred)<=abs(doy_event-doy_pred_2), doy_pred, doy_pred_2))
    rmse_photoperiod_allb=Metrics::rmse(pred_vs_actual_photoperiodb$doy_event, pred_vs_actual_photoperiodb$doy_pred_all)
    rmse_photoperiod_springb=Metrics::rmse(pred_vs_actual_photoperiodb$doy_event, pred_vs_actual_photoperiodb$doy_pred)
    rmse_photoperiod_fallb=Metrics::rmse(pred_vs_actual_photoperiodb$doy_event, pred_vs_actual_photoperiodb$doy_pred_2)
    
    ##photoperiodc
    pred_vs_actual_photoperiodc=clim_all_doy%>%
      rowwise()%>%
      filter(photoperiod>preds_photoperiodc&doy>0&doy<=173&!is.na(photoperiod))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), .groups='drop')%>%
      full_join(val_data%>%
                  select(!matches('preds')))%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 173, doy_pred))
    pred_vs_actual_photoperiod_2c=clim_all_doy%>%
      rowwise()%>%
      filter(photoperiod<preds_photoperiodc&doy>173&!is.na(photoperiod))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred_2=min(doy), .groups='drop')%>%
      full_join(val_data%>%
                  select(!matches('preds')))%>%
      mutate(doy_pred_2=ifelse(is.na(doy_pred_2), 173, doy_pred_2))
    pred_vs_actual_photoperiodc=full_join(pred_vs_actual_photoperiodc, pred_vs_actual_photoperiod_2c)%>%
      #pick the nearest one spring or fall
      mutate(doy_pred_all=ifelse(abs(doy_event-doy_pred)<=abs(doy_event-doy_pred_2), doy_pred, doy_pred_2))
    rmse_photoperiod_allc=Metrics::rmse(pred_vs_actual_photoperiodc$doy_event, pred_vs_actual_photoperiodc$doy_pred_all)
    rmse_photoperiod_springc=Metrics::rmse(pred_vs_actual_photoperiodc$doy_event, pred_vs_actual_photoperiodc$doy_pred)
    rmse_photoperiod_fallc=Metrics::rmse(pred_vs_actual_photoperiodc$doy_event, pred_vs_actual_photoperiodc$doy_pred_2)

    pred_vs_actual_MAT=clim_all_doy%>%
      select(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree, doy_event, preds_MAT)%>%
      distinct()%>%filter(!is.na(tree))
    rmse_MAT=Metrics::rmse(pred_vs_actual_MAT$doy_event, pred_vs_actual_MAT$preds_MAT)
    
    pred_vs_actual_MATa=clim_all_doy%>%
      select(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree, doy_event, preds_MATa)%>%
      distinct()%>%filter(!is.na(tree))
    rmse_MATa=Metrics::rmse(pred_vs_actual_MATa$doy_event, pred_vs_actual_MATa$preds_MATa)
    
    pred_vs_actual_lat=clim_all_doy%>%
      select(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree, doy_event, preds_lat)%>%
      distinct()%>%filter(!is.na(tree))
    rmse_lat=Metrics::rmse(pred_vs_actual_lat$doy_event, pred_vs_actual_lat$preds_lat)
    
    pred_vs_actual_lata=clim_all_doy%>%
      select(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree, doy_event, preds_lata)%>%
      distinct()%>%filter(!is.na(tree))
    rmse_lata=Metrics::rmse(pred_vs_actual_lata$doy_event, pred_vs_actual_lata$preds_lata)
    
    pred_vs_actual_forcea=clim_all_doy%>%
      rowwise()%>%
      filter(Force5s>preds_forcea&!is.na(Force5s))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), doy_event=unique(doy_event), .groups='drop')%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 365, doy_pred))
    rmse_forcea=Metrics::rmse(pred_vs_actual_forcea$doy_event, pred_vs_actual_forcea$doy_pred)
    
    pred_vs_actual_forceb=clim_all_doy%>%
      rowwise()%>%
      filter(Force5s>preds_forceb&!is.na(Force5s))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), doy_event=unique(doy_event), .groups='drop')%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 365, doy_pred))
    rmse_forceb=Metrics::rmse(pred_vs_actual_forceb$doy_event, pred_vs_actual_forceb$doy_pred)
    
    pred_vs_actual_forcec=clim_all_doy%>%
      rowwise()%>%
      filter(Force5s>preds_forcec&!is.na(Force5s))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), doy_event=unique(doy_event), .groups='drop')%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 365, doy_pred))
    rmse_forcec=Metrics::rmse(pred_vs_actual_forcec$doy_event, pred_vs_actual_forcec$doy_pred)
    
    pred_vs_actual_forced=clim_all_doy%>%
      rowwise()%>%
      filter(Force5s>preds_forced&!is.na(Force5s))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), doy_event=unique(doy_event), .groups='drop')%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 365, doy_pred))
    rmse_forced=Metrics::rmse(pred_vs_actual_forced$doy_event, pred_vs_actual_forced$doy_pred)
    
    pred_vs_actual_forcee=clim_all_doy%>%
      rowwise()%>%
      filter(Force5s>preds_forcee&!is.na(Force5s))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), doy_event=unique(doy_event), .groups='drop')%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 365, doy_pred))
    rmse_forcee=Metrics::rmse(pred_vs_actual_forcee$doy_event, pred_vs_actual_forcee$doy_pred)
    
    pred_vs_actual_forcef=clim_all_doy%>%
      rowwise()%>%
      filter(Force5s>preds_forcee&!is.na(Force5s))%>%
      group_by(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree)%>%
      summarise(doy_pred=min(doy), doy_event=unique(doy_event), .groups='drop')%>%
      mutate(doy_pred=ifelse(is.na(doy_pred), 365, doy_pred))
    rmse_forcef=Metrics::rmse(pred_vs_actual_forcef$doy_event, pred_vs_actual_forcef$doy_pred)
    
    pred_vs_actual_null=clim_all_doy%>%
      select(Species, Type, site, lat, lon, alt, Biome, year_of_event, tree, doy_event, preds_null)%>%
      distinct()%>%filter(!is.na(tree))
    rmse_null=Metrics::rmse(pred_vs_actual_null$doy_event, pred_vs_actual_null$preds_null)
    
    out[[k]]=list(sites2use=ifelse(by=='site', paste0(sites2use, collapse='|'), NA), 
                  years2use=ifelse(by=='year', paste0(years2use, collapse='|'), NA),
                  siteyrs2use=ifelse(by=='year_within_site',
                                   paste0(fit_data%>%select(site,year_of_event)%>%unite('siteyear', site, year_of_event)%>%
                                            pull(siteyear), collapse='|'), NA),
                  rmse_force=rmse_force, rmse_forcea=rmse_forcea, 
                  rmse_forceb=rmse_forceb, rmse_forcec=rmse_forcec, 
                  rmse_forced=rmse_forced,  rmse_forcee=rmse_forcee,
                  rmse_forcef=rmse_forcef,
                  rmse_chill=rmse_chill,
                  rmse_photoperiod_all=rmse_photoperiod_all, rmse_photoperiod_spring=rmse_photoperiod_spring,
                  rmse_photoperiod_fall=rmse_photoperiod_fall,
                  rmse_photoperiod_alla=rmse_photoperiod_alla, rmse_photoperiod_springa=rmse_photoperiod_springa,
                  rmse_photoperiod_falla=rmse_photoperiod_falla,
                  rmse_photoperiod_alla1=rmse_photoperiod_alla1, rmse_photoperiod_springa1=rmse_photoperiod_springa1,
                  rmse_photoperiod_falla1=rmse_photoperiod_falla1,
                  rmse_photoperiod_allb=rmse_photoperiod_allb, rmse_photoperiod_springb=rmse_photoperiod_springb,
                  rmse_photoperiod_fallb=rmse_photoperiod_fallb,
                  rmse_photoperiod_allc=rmse_photoperiod_allc, rmse_photoperiod_springc=rmse_photoperiod_springc,
                  rmse_photoperiod_fallc=rmse_photoperiod_fallc,
                  rmse_MAT=rmse_MAT, rmse_MATa=rmse_MATa, rmse_lat=rmse_lat,rmse_lata=rmse_lata,
                  rmse_null=rmse_null)
  }
  close(pb)
  return(out)
}

#run cross validations - takes a while
#lots of info messages on the join which you can ignore
# the pbar is somewhat helpful to know how long you have to go
# if using this more it would benefit from being parallelized...

cv_site=cv_forcing_chilling_photoperiod_CV(phen_data, clim_data, n_groups=100, prop_drop=0.1,
                                            by='site')
saveRDS (cv_site, 'calculations/cv_site.rds')


# 7 warning suggestsions on model identifiability
# could try scaling or just ignore...
#2 boundary fit warns, ignore for now
cv_year_within_site=cv_forcing_chilling_photoperiod_CV(phen_data, clim_data, n_groups=100,
                                                by='year_within_site', choose=FALSE)
saveRDS (cv_year_within_site, 'calculations/cv_year_within_site.rds')
#sce add progress bar for the impatient
#https://ryouready.wordpress.com/2009/03/16/r-monitor-function-progress-with-a-progress-bar/

#convert to df to plot
#there's surely a faster list indexing way to do this
#but this works
# it just make a ton of duplicates that you can then ditch with distinct()

#this for the density plots
cv_site_df_wide=data.table::rbindlist(cv_site, idcol='group')%>%
  select(-years2use, -siteyrs2use)%>%
  rename(null_val=rmse_null)%>%
  distinct()%>%  mutate_at(vars(contains('rmse')), .funs = rlang::list2(~. - null_val))%>%
  select(-null_val)%>%
  pivot_longer(cols=c(-group, -sites2use))%>%
  rename(model=name, rmse_delta=value)%>%
  mutate(model=gsub('rmse_', '', model))



#report only summary stats for chill as putting them in the figure
# they are so bad the x axis range has to be huge
mean(cv_site_df_wide$rmse_delta[cv_site_df_wide$model=='chill'])
quantile(cv_site_df_wide$rmse_delta[cv_site_df_wide$model=='chill'], c(0.05, 0.95))

# cv_year_df_wide=data.table::rbindlist(cv_year, idcol='group')%>%
#   select(-sites2use)%>%
#   rename(null_val=rmse_null)%>%
#   distinct()%>%
#   mutate_at(vars(contains('rmse')), .funs = rlang::list2(~. - null_val))%>%
#   select(-null_val)%>%
#   pivot_longer(cols=c(-group, -years2use))%>%
#   rename(model=name, rmse_delta=value)%>%
#   mutate(model=gsub('rmse_', '', model))



cv_year_within_site_df_wide=data.table::rbindlist(cv_year_within_site, idcol='group')%>%
  select(-sites2use, -years2use)%>%
  rename(null_val=rmse_null)%>%
  distinct()%>%
  mutate_at(vars(contains('rmse')), .funs = rlang::list2(~. - null_val))%>%
  select(-null_val)%>%
  pivot_longer(cols=c(-group, -siteyrs2use))%>%
  rename(model=name, rmse_delta=value)%>%
  mutate(model=gsub('rmse_', '', model))

#these so bad just report the summaries
mean(cv_year_within_site_df_wide$rmse_delta[cv_year_within_site_df_wide$model=='chill'])
quantile(cv_year_within_site_df_wide$rmse_delta[cv_year_within_site_df_wide$model=='chill'],
         c(0.05, 0.95))


#can look at percentage that are better or worse than null if desired
# cv_year_df_wide%>%
#   group_by(model)%>%
#   summarise(sum(rmse_delta>0)/dplyr::n())

View(cv_site_df_wide%>%
  group_by(model)%>%
  summarise(sum(rmse_delta>0)/dplyr::n()))

View(cv_year_within_site_df_wide%>%
  group_by(model)%>%
  summarise(prop_better_than_null=sum(rmse_delta>0)/dplyr::n()))



# #note the chill and fall photoperiod so bad you can't even get there
# null_year_model_plot = ggplot(data=cv_year_df_wide%>%
#                                 filter(!model%in%c('chill', 'photoperiod_fall',
#                                                    'force_lat', 'photoperiod_spring'))%>%
#                                 mutate(
#                                   model=case_when(
#                                     model=='force' ~ 'Forcing',
#                                     model=='lat' ~ 'Latitude',
#                                     model=='photoperiod_all' ~ 'Photoperiod',
#                                     model=='MAT' ~ 'MAT',
#                                     TRUE ~ 'NA')
#                                 ))+
#                       geom_density(aes(x=rmse_delta, group=model, fill=model, color=model), 
#                         alpha=0.5, adjust=2)+
#   theme(panel.background = element_blank())+
#   geom_vline(xintercept=0, linetype='dotted')+
#   facet_grid(rows = 'model')+
#   #xlim(-15,15)+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   labs(x = expression(delta[rmse(model - null)]))+
#   guides(fill=FALSE, color=FALSE)

#note the chill and fall photoperiod so bad you can't even get there

insets_year_all=cv_year_within_site_df_wide%>%
  group_by(model)%>%
  summarise(prop_worse_than_null=sum(rmse_delta>0)/dplyr::n())%>%
              filter(model%in%c('force', 'forcea',
                                'forceb', 'forcee',
                                'photoperiod_all',
                                'photoperiod_alla',
                                'photoperiod_allb',
                                'lat', 'MAT'))%>%
              mutate(
                model=case_when(
                  model=='force' ~ 'forcing',
                  model=='forcea' ~ 'forcing by lat',
                  model=='forceb' ~ 'forcing by chill',
                  model=='forcee' ~ 'forcing by chill and lat',
                  model=='lat' ~ 'latitude',
                  model=='photoperiod_all' ~ 'photoperiod',
                  model=='photoperiod_alla' ~ 'photoperiod by MAT',
                  model=='photoperiod_allb' ~ 'photoperiod by lat',
                  model=='MAT' ~ 'MAT',
                  TRUE ~ 'NA')
              )

null_year_model_plot = ggplot(data=cv_year_within_site_df_wide%>%
                                filter(model%in%c('force', 'forcea',
                                                  'forceb', 'forcee',
                                                  'photoperiod_all',
                                                  'photoperiod_alla',
                                                  'photoperiod_allb',
                                                  'lat', 'MAT'))%>%
                                mutate(
                                  model=case_when(
                                    model=='force' ~ 'forcing',
                                    model=='forcea' ~ 'forcing by lat',
                                    model=='forceb' ~ 'forcing by chill',
                                    model=='forcee' ~ 'forcing by chill and lat',
                                    model=='lat' ~ 'latitude',
                                    model=='photoperiod_all' ~ 'photoperiod',
                                    model=='photoperiod_alla' ~ 'photoperiod by MAT',
                                    model=='photoperiod_allb' ~ 'photoperiod by lat',
                                    model=='MAT' ~ 'MAT',
                                    TRUE ~ 'NA')
                                )
                               )+
  geom_density(aes(x=rmse_delta, group=model, fill=model, color=model), 
               alpha=0.5, adjust=2)+
  theme(panel.background = element_blank())+
  geom_vline(xintercept=0, linetype='dotted')+
  #facet_grid(rows = 'model')+
  facet_wrap(~model, ncol = 2)+
  xlim(-15,15)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = expression(delta[rmse(model - null)]))+
  guides(fill=FALSE, color=FALSE)+
  geom_text(data=insets, aes(x=10, y=4, label=prop_worse_than_null), size=3)

ggsave("plots/null_year_model_plot.jpg", null_year_model_plot2, height=12, width=6)

##
insets_year_subset=cv_year_within_site_df_wide%>%
  group_by(model)%>%
  summarise(prop_worse_than_null=sum(rmse_delta>0)/dplyr::n())%>%
  filter(model%in%c('force', 
                      'photoperiod_all',
                      'lat', 'MAT'))%>%
  mutate(
    model=case_when(
      model=='force' ~ 'forcing',
      model=='lat' ~ 'latitude',
      model=='photoperiod_all' ~ 'photoperiod',
      model=='MAT' ~ 'MAT',
      TRUE ~ 'NA')
  )

#just select subset for paper
null_year_model_plot_subset = ggplot(data=cv_year_within_site_df_wide%>%
                                 filter(model%in%c('force', 
                                                   'photoperiod_all',
                                                   'lat', 'MAT'))%>%
                                 mutate(
                                   model=case_when(
                                     model=='force' ~ 'forcing',
                                     model=='lat' ~ 'latitude',
                                     model=='photoperiod_all' ~ 'photoperiod',
                                     model=='MAT' ~ 'MAT',
                                     TRUE ~ 'NA')
                                 )
)+
  geom_density(aes(x=rmse_delta, group=model, fill=model, color=model), 
               alpha=0.5, adjust=2)+
  theme(panel.background = element_blank())+
  geom_vline(xintercept=0, linetype='dotted')+
  facet_grid(rows = 'model')+
  #facet_wrap(~model, ncol = 1)+
  xlim(-15,15)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = expression(delta[rmse(model - null)]))+
  guides(fill=FALSE, color=FALSE)+
  geom_text(data=insets_year_subset, aes(x=10, y=4, label=prop_worse_than_null), size=3)

ggsave("plots/null_year_model_plot_subset.jpg",null_year_model_plot_subset, height=12, width=6)

insets_site_all=cv_site_df_wide%>%
  group_by(model)%>%
  summarise(prop_worse_than_null=sum(rmse_delta>0)/dplyr::n())%>%
  filter(model%in%c('force', 'forcea',
                    'forceb', 'forcee',
                    'photoperiod_all',
                    'photoperiod_alla',
                    'photoperiod_allb',
                    'lat', 'MAT'))%>%
  mutate(
    model=case_when(
      model=='force' ~ 'forcing',
      model=='forcea' ~ 'forcing by lat',
      model=='forceb' ~ 'forcing by chill',
      model=='forcee' ~ 'forcing by chill and lat',
      model=='lat' ~ 'latitude',
      model=='photoperiod_all' ~ 'photoperiod',
      model=='photoperiod_alla' ~ 'photoperiod by MAT',
      model=='photoperiod_allb' ~ 'photoperiod by lat',
      model=='MAT' ~ 'MAT',
      TRUE ~ 'NA')
  )

null_site_model_plot = ggplot(data=cv_site_df_wide%>%
  filter(model%in%c('force', 'forcea',
                    'forceb', 'forcee',
                    'photoperiod_all',
                    'photoperiod_alla',
                    'photoperiod_allb',
                    'lat', 'MAT'))%>%
  mutate(
    model=case_when(
      model=='force' ~ 'forcing',
      model=='forcea' ~ 'forcing by lat',
      model=='forceb' ~ 'forcing by chill',
      model=='forcee' ~ 'forcing by chill and lat',
      model=='lat' ~ 'latitude',
      model=='photoperiod_all' ~ 'photoperiod',
      model=='photoperiod_alla' ~ 'photoperiod by MAT',
      model=='photoperiod_allb' ~ 'photoperiod by lat',
      model=='MAT' ~ 'MAT',
      TRUE ~ 'NA')))+
  geom_density(aes(x=rmse_delta, group=model, fill=model, color=model), 
               alpha=0.5, adjust=2)+
  theme(panel.background = element_blank())+
  geom_vline(xintercept=0, linetype='dotted')+
  facet_grid(rows = 'model')+
  #facet_wrap(~model, ncol = 2)+
  xlim(-20,20)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = expression(delta[rmse(model - null)]))+
  guides(fill=FALSE, color=FALSE)+
  geom_text(data=insets_site_all, aes(x=10, y=0.3, label=prop_worse_than_null), size=3)


insets_site_subset=cv_site_df_wide%>%
  group_by(model)%>%
  summarise(prop_worse_than_null=sum(rmse_delta>0)/dplyr::n())%>%
  filter(model%in%c('force', 
                    'photoperiod_all',
                    'lat', 'MAT'))%>%
  mutate(
    model=case_when(
      model=='force' ~ 'forcing',
      model=='lat' ~ 'latitude',
      model=='photoperiod_all' ~ 'photoperiod',
      model=='MAT' ~ 'MAT',
      TRUE ~ 'NA')
  )

#just select subset for paper
null_site_model_plot_subset = ggplot(data=cv_site_df_wide%>%
                                       filter(model%in%c('force', 
                                                         'photoperiod_all',
                                                         'lat', 'MAT'))%>%
                                       mutate(
                                         model=case_when(
                                           model=='force' ~ 'forcing',
                                           model=='lat' ~ 'latitude',
                                           model=='photoperiod_all' ~ 'photoperiod',
                                           model=='MAT' ~ 'MAT',
                                           TRUE ~ 'NA')
                                       )
)+
  geom_density(aes(x=rmse_delta, group=model, fill=model, color=model), 
               alpha=0.5, adjust=2)+
  theme(panel.background = element_blank())+
  geom_vline(xintercept=0, linetype='dotted')+
  facet_grid(rows = 'model')+
  #facet_wrap(model, ncol = 1)+
  xlim(-20,20)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = expression(delta[rmse(model - null)]))+
  guides(fill=FALSE, color=FALSE)+
  geom_text(data=insets_site_subset, aes(x=10, y=0.3, label=prop_worse_than_null), size=3)


#combine plots
joint_plot=
  egg::ggarrange(null_site_model_plot_subset+
                   ggtitle("spatial")+
                   theme(plot.title = element_text(hjust = 0.5),
                         axis.title.x = element_blank(),
                         strip.background = element_blank(),
                         strip.text.y = element_blank()),
                 null_year_model_plot_subset + 
                   ggtitle("temporal")+
                   theme(plot.title = element_text(hjust = 0.5))+
                   theme(#axis.text.y = element_blank(),
                         #axis.ticks.y = element_blank(),
                         axis.title.y = element_blank(),
                         axis.title.x = element_blank(),
                         strip.text.y = element_text(size = 14)),
                 nrow = 1, labels = c('a', 'b'),
bottom = grid::textGrob(
  #this works
  #paste0('delta \n', expression(better %<->% worse)),
  #this does not work
  #paste0(expression(delta[rmse(model - null)]), '\n', expression(better %<->% worse)),
  expression(atop(delta[rmse(model - null)], better %<->% worse)),
  #expression(better %<->% worse),
  #expression(delta[rmse(model - null)] \n better %<->% worse),
  #expression(paste('delta \n', better %<->% worse)),
  gp = grid::gpar(fontface = 3, fontsize = 18),
  #hjust = 1,
  just='center'
  #x = 0.1
))

ggsave("plots/joint_plot.jpg", joint_plot, height=12, width=8)




#add intuitive label axis
spatial_cv=gridExtra::grid.arrange(
  null_site_model_plot,
  nrow = 1,
  top = "Predictive performace (spatial variation)",
  bottom = textGrob(
    paste0('delta \n', expression(better %<->% worse)),
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

