
# library(brms) # for the analysis
# library(rstan)
# library(tidybayes)
library(sjstats)
# library(sjmisc)
# library(sjPlot)
library(tidyverse)
library(broom)
library(magrittr)
library(taxonlookup)
library(rjson)
library(sp)
library(GSIF)
library(AICcmodavg)
library(lmerTest)
library(medfate)
# Rcpp::sourceCpp('soil.cpp')
# library(seewave)
library(purrr)
library(future)
library(furrr)
plan('multicore')
options('future.global.maxsize'=2*1024*1024^2)

## 0.1 Plant names to be modelized ---------------

path <- "~/Chapter_3/data/species_daylight/"
species_names <- list.files(path = path)

purrr::map(species_names,function(.x){
  
  
  ## 1.0 Data load -----------------------
  
  load(paste0(path,.x))
  
  print(.x)
  
  x2 %>% select(pl_code,pl_dbh,pl_sapw_area,pl_species,si_code) %>% unique() -> metadata
  
  metadata <- metadata %>% mutate(pl_species = fct_recode(pl_species, 
                                                          "Vouacapoua americana"="Vacapoua americana"),
                                  pl_species = as.character(pl_species))
  
  metadata %>%
    pull(pl_species) %>%
    unique() %>%
    lookup_table(missing_action = 'NA', by_species = TRUE) %>%
    rownames_to_column('pl_species') %>%
    left_join(metadata, ., by = 'pl_species') -> res
  
  return(res)
  
}) %>% bind_rows() -> data_dbh_sapwood

data_dbh_sapwood %>% mutate(ba = (pl_dbh/2)^2*pi,
                            pl_dbh=log(pl_dbh),
                            pl_sapw_area = log(pl_sapw_area),
                            ba = log(ba)) %>%
  filter(pl_sapw_area<(-0.2)) %>% select(pl_code) %>% unique() -> filter_sites

save(filter_sites, file = 'filter_sites_sapw_units_wrong.Rdata')

data_dbh_sapwood %>%  
  mutate(pl_sapw_area = case_when(pl_code %in% filter_sites$pl_code ~ pl_sapw_area*10000,
                            !pl_code %in% filter_sites$pl_code ~ pl_sapw_area),
         ba = (pl_dbh/2)^2*pi,
         pl_dbh=log(pl_dbh),
         log_pl_sapw_area = log(pl_sapw_area),
         log_ba = log(ba)) -> data_dbh_sapwood
  # ggplot(aes(log_ba,log_pl_sapw_area,color=group))+
  # geom_point()

lm(log_pl_sapw_area~log_ba*group,data_dbh_sapwood)->mod
summary(mod)

saveRDS(mod,file="dbh_sapw_area_model.rds")
