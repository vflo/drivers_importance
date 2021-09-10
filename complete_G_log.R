library(tidyverse)
library(taxonlookup)
library(broom)
library(magrittr)
library(rjson)
library(sp)
library(lme4)
source("~/Drivers_importance/soil_functions.R")
library(purrr)
library(future)
library(furrr)
library(partR2)
plan('multisession')
options('future.global.maxsize'=2*1024*1024^2)

num_NA <- function(x, ...) {
  if (length(x) == 0)
    return(NA)
  return(x)
}


## 0.1 Plant names to be modelized ---------------

path <- "~/Drivers_importance/data/site_daylight/"
site_names <- list.files(path = path)

source('~/Drivers_importance/treatment_filter_only_control.R')

load("~/Drivers_importance/data/swc_ERA5_land.RData") #nearest cell extract ERA5 land 9x9km 1980 to present

ERA5_land_data <- swc_ERA5_land %>% 
  rename(ERA5_swc = swc_ERA5_land,
         TIMESTAMP = TIMESTAMP_daylight) %>% 
  mutate(TIMESTAMP = lubridate::as_date(TIMESTAMP))

dbh_sapw_mod <- readRDS(file="~/Drivers_importance/dbh_sapw_area_model.rds")

## import PET
PET <- read_csv('~/Drivers_importance/data/PET.csv') %>% 
  dplyr::rename(si_lat = lati,si_long = long)
BIOS <- read_csv("~/Drivers_importance/data/CHELSA_BIOS.csv") %>% mutate(si_long = lon, si_lat = lat)

load("~/Drivers_importance/data/sand.RData")
load("~/Drivers_importance/data/clay.RData")
load("~/Drivers_importance/data/elevation.RData")

load("~/Drivers_importance/data/sw_ERA5_18.RData")
load("~/Drivers_importance/data/sw_ERA5_6.RData")


sw <- sw_ERA5_18 %>% 
  left_join(sw_ERA5_6) %>% 
  mutate(sw = sw_ERA5_18 - sw_ERA5_6,
         si_code = as.character(si_code)) %>% 
  rename("TIMESTAMP" = "TIMESTAMP_daylight")%>% 
  mutate(TIMESTAMP = lubridate::as_date(TIMESTAMP))

rm(sw_ERA5_18)
rm(sw_ERA5_6)


site_names <- sample(site_names)


## 0.2 MODELIZATION ----------------------

furrr::future_map(site_names,.progress=TRUE, .f=function(.x){

    ## 1.0 Data load -----------------------
  
  load(paste0(path,.x))
  
  print(.x)
  
  ## 1.1 Conditionals to allow model calculation -----
  
  x2 <- x2 %>% mutate(st_treatment = case_when(is.na(st_treatment) ~ "None",
                                               !is.na(st_treatment) ~ st_treatment))
  if(any(colnames(x2) == "vpd_mean")
  ){ 
    
    if(!any(colnames(x2) == "ppfd_in_mean")){x2 <- x2 %>% mutate(ppfd_in_mean = c(NA))}
    if(!any(colnames(x2) == "sw_in_mean")){x2 <- x2 %>% mutate(sw_in_mean = c(NA))}
    if(!any(colnames(x2) == "swc_shallow_mean")){x2 <- x2 %>% mutate(swc_shallow_mean = c(NA))}
    
    metadata <- x2 %>%
      mutate(pl_species = as.character(pl_species),
             TIMESTAMP = lubridate::ymd(TIMESTAMP))
    
    metadata %>%
      pull(pl_species) %>%
      unique() %>%
      lookup_table(missing_action = 'NA', by_species = TRUE) %>%
      rownames_to_column('pl_species') %>%
      left_join(metadata, ., by = 'pl_species') -> x2
    

    x2 <- x2 %>% cbind(elev = elevation %>% 
                         filter(si_code == unique(x2$si_code)) %>% 
                         pull(elevation))
    
    
    ## 2.1 Data preparation ---------------
    
    x2 %>% left_join(ERA5_land_data, by = c("TIMESTAMP", "si_code")) %>%
      mutate(swc = ERA5_swc) %>% 
      filter(!is.na(vpd_mean)) %>% 
      left_join(sw,by = c("TIMESTAMP", "si_code")) %>% 
      mutate(TIMESTAMP = lubridate::ymd(TIMESTAMP))%>% 
      mutate(
        ba = ((pl_dbh/2)^2)*pi,
        log_ba = log(ba),
        log_pl_sapw_area = predict.lm(object = dbh_sapw_mod, tibble(log_ba,group)),
        pl_sapw_area = case_when(is.na(pl_sapw_area) ~ exp(log_pl_sapw_area),
                                 !is.na(pl_sapw_area) ~ pl_sapw_area),
        sapflow_mean = case_when(pl_sap_units == "cm3 cm-2 h-1" ~ sapflow_mean,# cm3 cm-2 h-1
                                 pl_sap_units == "cm3 h-1" ~ sapflow_mean/pl_sapw_area),#from cm3 h-1 to cm3 cm2 h-1
        si_elev = coalesce(si_elev,elev),
        ppfd_ERA = LakeMetabolizer::sw.to.par.base(sw),
        ppfd_sw_in = LakeMetabolizer::sw.to.par.base(sw_in_mean),
        ppfd = case_when(
          # We priorice ppfd calculated with on-site sw, and secondly ERA5 sw 
          is.na(sw_in_mean)~ppfd_ERA, 
          !is.na(sw_in_mean)~as.double(ppfd_sw_in),
          !is.na(ppfd_in_mean)~as.double(ppfd_in_mean)
          )
      ) -> x2
    
    
    x2%>% 
      filter(st_treatment %in% treatments,
             !is.na(vpd_mean),
             !is.na(sapflow_mean),
             !is.na(ta_mean),
      )%>%
      mutate(diff = swc / lag(swc, default = first(swc), order_by = TIMESTAMP)) %>% #slope to obtain swc recharge
      filter(#lag(diff,n=1)<=1,
             #lead(diff,n=1)<=1,
             diff<=1
      ) %>% #filter no rain days + - 1 day
      dplyr::mutate(
        sapflow_mean = case_when(pl_sens_meth == "HD" & !pl_sens_calib ~ 
                                   sapflow_mean+0.40488*sapflow_mean,# calibration correction from flo et al. 2019
                                 pl_sens_meth == "HD" & is.na(pl_sens_calib) ~ 
                                   sapflow_mean+0.40488*sapflow_mean,# calibration correction from flo et al. 2019
                                 pl_sens_meth == "HD" & pl_sens_calib ~ sapflow_mean,
                                 pl_sens_meth != "HD" ~ sapflow_mean),
        cond_coeff = 1.158e2 + 4.236e-1 * ta_mean, #[kPa m3 kg-1]
        conductance = cond_coeff*(sapflow_mean/(60*60*10^3))/vpd_mean,#cm3 cm-2 h-1 to m3 vapor cm-2 s-1
        conductance = case_when(!is.na(si_elev) ~ conductance*44.6*(273/(273+ta_mean))*
                                  (101.325*exp(-0.00012*si_elev)/101.325),
                                is.na(si_elev) ~ conductance*44.6*(273/(273+ta_mean))),#from m3 cm-2 s-1 to mol*cm-2 s-1
        G_sw = conductance*10^4,#from mol cm-2 s-1 to mol m-2 s-1
        G_max = quantile(conductance, probs = 0.95,na.rm = TRUE),
        G_min = min(conductance,na.rm = TRUE),
        G_range = G_max - G_min,
        hour = lubridate::hour(TIMESTAMP),
        DOY = lubridate::round_date(TIMESTAMP,unit = 'day') %>% as.numeric(),
        si_lat = mean(si_lat, na.rm =TRUE),
        si_long = mean(si_long, na.rm = TRUE),
        st_sand_perc = mean(st_sand_perc, na.rm = TRUE),
        st_clay_perc = mean(st_clay_perc, na.rm = TRUE),
        st_silt_perc = mean(st_silt_perc, na.rm = TRUE),
        pl_height = coalesce(pl_height, st_height)
      ) %>%
      arrange(TIMESTAMP) -> faa
    
    
  }else{faa <- tibble(swc = NA)}
  
  if (!is.null(faa) & any(!is.na(faa$swc))){
    
    # 2.2.1 Soil data import --------------

    faa <- faa %>% cbind(sand = sand %>% 
                           filter(si_code == unique(faa$si_code)) %>% 
                           pull(sand))
    faa <- faa %>% cbind(clay = clay %>% 
                           filter(si_code == unique(faa$si_code)) %>% 
                           pull(clay))
    
    ## 2.2.2 Soil WWP and FC calculation ---------
    
    swp <- NA
    swp.sat <- NA
    if(any(is.na(unique(faa$st_sand_perc)),is.na(unique(faa$st_clay_perc)))){
      swp <- clapp.and.hornberger(faa$swc,
                                  percent.clay=unique(faa$clay),
                                  percent.sand=unique(faa$sand))[["suction"]] #mm
      swp.sat <- clapp.and.hornberger(faa$swc,
                                      percent.clay=unique(faa$clay),
                                      percent.sand=unique(faa$sand))[["suction.sat"]] #mm
      
      swp <- -(swp/101972) #MPa
      swp.sat <- -(swp.sat/101972) #MPa
    }
    if(all(!is.na(unique(faa$st_sand_perc)),!is.na(unique(faa$st_clay_perc)))){
      swp <- clapp.and.hornberger(faa$swc,
                                  percent.clay=unique(faa$st_clay_perc),
                                  percent.sand=unique(faa$st_sand_perc))[["suction"]] #mm
      swp.sat <- clapp.and.hornberger(faa$swc,
                                      percent.clay=unique(faa$st_clay_perc),
                                      percent.sand=unique(faa$st_sand_perc))[["suction.sat"]] #mm
      swp <- -(swp/101972) #MPa
      swp.sat <- -(swp.sat/101972) #MPa
    }
    
    faa <- faa %>% cbind(data.frame(swp = swp,swp.sat = swp.sat), row.names = NULL)
    
    save(faa, file = paste0('./data/models_raw/complete_G_log/',.x))
    
    
    ## 2.3 Environmental filters ---------
    faa%>%
      filter(!is.na(G_sw),
             !is.na(swc),
             !is.na(vpd_mean),
             !is.na(ppfd),
             G_sw > 0,
             vpd_mean >= 0.3,
             swc>0,
             ppfd>0,
             swc <= 0.3
      )-> faa
    
    if (!is.null(faa) & any(!is.na(faa$swc))){
    
    #Metadata and biogeo variables
    faa %>%
      left_join(PET %>% dplyr::select(-c(si_lat,si_long)),by=c('si_code'))%>%
      left_join(BIOS %>% dplyr::select(-c(si_lat,si_long)),by=c('si_code'))%>%
      group_by(pl_code) %>% 
      mutate(pl_height = coalesce(pl_height, st_height),
             n = n()) %>%
      ungroup() %>% 
      summarise(si_code = unique(si_code),
                pl_dbh = weighted.mean(pl_dbh,n),
                pl_height = weighted.mean(pl_height,n),
                si_long = unique(si_long),
                si_lat = unique(si_lat),
                si_elev = unique(si_elev),
                MAP = unique(MAP),
                MAT = unique(MAT)/10,
                PET = unique(PET),
                PPET = MAP/PET,
                si_biome = unique(si_biome),
                min_timestamp = min(TIMESTAMP,na.rm=TRUE),
                max_timestamp = max(TIMESTAMP,na.rm=TRUE),
                pl_sap_units = unique(pl_sap_units),
                swc_shallow = ifelse(all(is.na(swc_shallow_mean)),
                                     NA_integer_,
                                     mean(swc_shallow_mean,na.rm = TRUE)
                                     ),
                sw_in_mean = ifelse(all(is.na(sw_in_mean)),
                                   NA_integer_,
                                   mean(sw_in_mean,na.rm = TRUE)
                                   ),
                ppfd_in_mean = ifelse(all(is.na(ppfd_in_mean)),
                                   NA_integer_,
                                   mean(ppfd_in_mean,na.rm = TRUE)
                                   ),
                ppfd_sw_in = ifelse(all(is.na(ppfd_sw_in)),
                                      NA_integer_,
                                      mean(ppfd_sw_in,na.rm = TRUE)
                                    ),
                ppfd = ifelse(all(is.na(ppfd)),
                                    NA_integer_,
                                    mean(ppfd,na.rm = TRUE)
                              )
      )-> env_data
    
    # save(faa, file = paste0('./data/models_raw/complete_G_log/',.x))
    
    ## 3.0 models -----------------  
    
    #plant level filtering
    faa %>%
        group_by(pl_code) %>% 
        mutate(log_vpd_mean = -log(vpd_mean),
               log_swc = log(swc),
               log_ppfd = log(ppfd/1000),
               min_swc = min(swc,na.rm = TRUE),
               max_swc = max(swc,na.rm = TRUE),
               range_swc = max_swc - min_swc,
               min_vpd = min(vpd_mean,na.rm = TRUE),
               max_vpd = max(vpd_mean,na.rm = TRUE),
               range_vpd = max_vpd - min_vpd,
               min_ppfd = min(ppfd,na.rm = TRUE),
               max_ppfd = max(ppfd,na.rm = TRUE),
               range_ppfd = max_ppfd - min_ppfd,
               n_day_plant = n_distinct(vpd_mean))%>%
        filter(n_day_plant>=10,
               range_swc>0.05|range_vpd>0.5) %>% 
        ungroup() %>%
        mutate(min_day_plant = min(n_day_plant)) ->faa2
    
    #site level minimum data conditional
    if (all(n_distinct(faa2$vpd_mean) >= 9,
          n_distinct(faa2$swc) >= 9,
          n_distinct(faa2$ppfd) >= 9,
          max(faa2$max_ppfd) >= 400)){

        test_1 <- unique(faa2$pl_species) %>% length()>1
        test_2 <- faa2 %>% 
          group_by(pl_species) %>% 
          summarise(n_dis = n_distinct(pl_code)) %>%
          mutate(logi=any(n_dis > 1)) %>%
          dplyr::select(logi) %>% 
          unique() %>% 
          as.logical()

        
        if(test_1 & test_2){

          type = 1

            model_G_sw_log <- try(lmer(G_sw~log_vpd_mean+log_swc+log_ppfd+
                                         (-1+log_vpd_mean+log_swc+log_ppfd|pl_species)+
                                         (1|pl_species:pl_code),REML = TRUE,
                                       data = faa2))
            model_G_sw_log_vpd <- try(lmer(G_sw~log_vpd_mean+
                                             (log_vpd_mean|pl_species)+
                                             (1|pl_species:pl_code),REML = TRUE,data = faa2))
            model_G_sw_log_swc <- try(lmer(G_sw~log_swc+
                                             (log_swc|pl_species)+
                                             (1|pl_species:pl_code),REML = TRUE,data = faa2))
            model_G_sw_log_ppfd <- try(lmer(G_sw~log_ppfd+
                                              (log_ppfd|pl_species)+
                                              (1|pl_species:pl_code),REML = TRUE,data = faa2))
            relgrad <- with(model_G_sw_log@optinfo$derivs,solve(Hessian,gradient))

            if(max(abs(relgrad))>0.0001){
              test_1 <- FALSE
              test_2 <- TRUE}

      }

        
        if(test_1 & !test_2){
          
          type = 2

          model_G_sw_log <- try(lmer(G_sw~log_vpd_mean+log_swc+log_ppfd+
                                       (log_vpd_mean+log_swc+log_ppfd|pl_species),REML = TRUE,
                                     data = faa2))
          model_G_sw_log_vpd <- try(lmer(G_sw~log_vpd_mean+
                                       (log_vpd_mean|pl_species),REML = TRUE,
                                       data = faa2))
          model_G_sw_log_swc <- try(lmer(G_sw~log_swc+
                                       (log_swc|pl_species),REML = TRUE,
                                       data = faa2))
          model_G_sw_log_ppfd <- try(lmer(G_sw~log_ppfd+
                                       (log_ppfd|pl_species),REML = TRUE,
                                       data = faa2))

        }
        
        if(!test_1 & test_2){

           type = 3

             model_G_sw_log <- try(lmer(G_sw~log_vpd_mean+log_swc+log_ppfd+
                                          (1|pl_code),REML = TRUE,
                                        data = faa2))
             model_G_sw_log_vpd <- try(lmer(G_sw~log_vpd_mean+
                                          (1|pl_code),REML = TRUE,
                                        data = faa2))
             model_G_sw_log_swc <- try(lmer(G_sw~log_swc+
                                          (1|pl_code),REML = TRUE,
                                        data = faa2))
             model_G_sw_log_ppfd <- try(lmer(G_sw~log_ppfd+
                                          (1|pl_code),REML = TRUE,
                                        data = faa2))
             }  
  
        if(!test_1 & !test_2){

          type = 4

          model_G_sw_log <- try(lm(G_sw~log_vpd_mean+log_swc+log_ppfd,data = faa2))
          model_G_sw_log_vpd <- try(lm(G_sw~log_vpd_mean,data = faa2))
          model_G_sw_log_swc <- try(lm(G_sw~log_swc,data = faa2))
          model_G_sw_log_ppfd <- try(lm(G_sw~log_ppfd,data = faa2))

        } 

        # partR2(model_G_sw_log,
        #        partvars = c("log_vpd_mean", "log_swc", "log_ppfd"),
        #        R2_type = "conditional", nboot = 100, CI = 0.95,
        #        data = faa2)
        # r2glmm::r2beta(model_G_sw_log,method="nsj")

          r2_G_log <- ifelse(class(model_G_sw_log) == "try-error", as.numeric("NA"),
                             model_G_sw_log %>% MuMIn::r.squaredGLMM() %>% .[1,"R2c"])
          
          r2_G_log_mar <- ifelse(class(model_G_sw_log) == "try-error", as.numeric("NA"),
                             model_G_sw_log %>% MuMIn::r.squaredGLMM() %>% .[1,"R2m"])

          rel_r2 <- r2glmm::r2beta(model_G_sw_log, method="nsj")
          rel_r2 <- rel_r2 %>% as_tibble() %>% dplyr::select(Effect,Rsq)
          vpd_rel = rel_r2[which(rel_r2$Effect=='log_vpd_mean'),2]$Rsq
          swc_rel = rel_r2[which(rel_r2$Effect=='log_swc'),2]$Rsq
          ppfd_rel = rel_r2[which(rel_r2$Effect=='log_ppfd'),2]$Rsq

          relat = vpd_rel+swc_rel+ppfd_rel
          vpd_rel <- vpd_rel/relat
          swc_rel <- swc_rel/relat
          ppfd_rel <- ppfd_rel/relat

          r2_G_log_vpd <- ifelse(class(model_G_sw_log_vpd) == "try-error", as.numeric("NA"),
                                 model_G_sw_log_vpd %>% MuMIn::r.squaredGLMM() %>% .[1,"R2c"])

          r2_G_log_swc <- ifelse(class(model_G_sw_log_swc) == "try-error", as.numeric("NA"),
                                 model_G_sw_log_swc %>% MuMIn::r.squaredGLMM() %>% .[1,"R2c"])

          r2_G_log_ppfd <- ifelse(class(model_G_sw_log_ppfd) == "try-error", as.numeric("NA"),
                                  model_G_sw_log_ppfd %>% MuMIn::r.squaredGLMM() %>% .[1,"R2c"])
          
          r2_G_log_vpd_mar <- ifelse(class(model_G_sw_log_vpd) == "try-error", as.numeric("NA"),
                                 model_G_sw_log_vpd %>% MuMIn::r.squaredGLMM() %>% .[1,"R2m"])
          
          r2_G_log_swc_mar <- ifelse(class(model_G_sw_log_swc) == "try-error", as.numeric("NA"),
                                 model_G_sw_log_swc %>% MuMIn::r.squaredGLMM() %>% .[1,"R2m"])
          
          r2_G_log_ppfd_mar <- ifelse(class(model_G_sw_log_ppfd) == "try-error", as.numeric("NA"),
                                  model_G_sw_log_ppfd %>% MuMIn::r.squaredGLMM() %>% .[1,"R2m"])

          df <- tibble(n_days_complete = nrow(faa2),
                       r2_G_log = r2_G_log,
                       r2_G_log_vpd = r2_G_log_vpd,
                       r2_G_log_swc = r2_G_log_swc,
                       r2_G_log_ppfd = r2_G_log_ppfd,
                       r2_G_log_mar = r2_G_log_mar,
                       r2_G_log_vpd_mar = r2_G_log_vpd_mar,
                       r2_G_log_swc_mar = r2_G_log_swc_mar,
                       r2_G_log_ppfd_mar = r2_G_log_ppfd_mar,
                       model_type_complete = type,
                       relat_sum = relat,
                       vpd_rel = vpd_rel,
                       swc_rel = swc_rel,
                       ppfd_rel = ppfd_rel)
# 
#         return(df)
#           
#         }) %>% bind_rows() -> df_res

        
        print(unique(x2$si_code))
        
        model <- list(df, 
                      env_data,
                      faa2,
                      model_G_sw_log
                      )
        
        save(model,file=paste0(getwd(),"/data/models/complete_G_log/",unique(x2$si_code),".RData"))
      }else{NULL}
      
    }else{NULL}
    
  }else{NULL}
  
  rm(x2)
  gc()
  
}
)
