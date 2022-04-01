library(tidyverse)
library(lmerTest)
library(cowplot)
library(purrr)
library(future)
library(furrr)
plan('multisession')
options('future.global.maxsize'=2*1024*1024^2)

custom_palette <- c("#e6545d", "#f49f27", "#ACE4AA", "#2b99ca", "#40596f")

## 0.1 Site names  ---------------

path <- "~/Drivers_importance/data/site_daylight/"
site_names <- list.files(path = path)

source('~/Drivers_importance/treatment_filter_only_control.R')
biomes <- read.csv("data/biome.csv")
load("~/Drivers_importance/data/swc_ERA5_land.RData") #nearest cell extract ERA5 land 9x9km 1980 to present

ERA5_land_data <- swc_ERA5_land %>% 
  rename(ERA5_swc = swc_ERA5_land,
         TIMESTAMP = TIMESTAMP_daylight) %>% 
  mutate(TIMESTAMP = lubridate::as_date(TIMESTAMP))


site_names <- sample(site_names)



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
    
    ## Join ERA5 swc
    
    x2 %>% left_join(ERA5_land_data, by = c("TIMESTAMP", "si_code")) %>%
      left_join(biomes, by = "si_code") %>% 
      mutate(swc = ERA5_swc) %>% 
      dplyr::select(TIMESTAMP, si_code, swc, swc_shallow_mean,biome) %>% 
      distinct()
    
  }
}
) %>% bind_rows()-> swc_df

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

swc_mod <- lmer(swc_shallow_mean~ swc+ (swc|si_code), data = swc_df)
summary(swc_mod)
# swc_mod <- lmer(swc~swc_shallow_mean + (1|si_code), data = swc_df)
# summary(swc_mod)
r_sqr <- MuMIn::r.squaredGLMM(swc_mod)

predict_swc <- predict(swc_mod, re.form =(~swc|si_code))


(swc_df %>%  
  na.omit() %>% cbind(predict_swc) %>%
  mutate(density = get_density(swc_shallow_mean,predict_swc, n = 200)) %>%
  ggplot(aes(swc_shallow_mean, predict_swc,color =density))+
  geom_point(alpha=0.5)+
  geom_abline(slope = 1, intercept=0)+
  # geom_abline(slope=swc_mod@beta[2], intercept = swc_mod@beta[1],
  #             color="blue", linetype="dashed", size=1)+
  scale_color_gradient(low = "#40596f", high = "#ACE4AA", na.value = NA)+
  # viridis::scale_color_viridis()+
  labs(color = "Density")+
  annotate("text", x = 0.55, y = 0.1, 
           label = "paste(italic(R) ^ 2, conditional,\" = 0.91\")", parse = TRUE)+
  theme_bw()+
  ylab("SWC predicted")+
  xlab("SAPFLUXNET SWC") -> fig1)


# swc_df %>% 
#   na.omit() %>%
#   group_by(si_code) %>% 
#   summarise("Pearson's correlation" = cor(swc,swc_shallow_mean)) %>% 
#   ggplot(aes(x=`Pearson's correlation`)) + 
#   geom_histogram(color="black", fill="white")+ 
#   geom_vline(aes(xintercept=median(`Pearson's correlation`)),
#              color="blue", linetype="dashed", size=1)+
#   theme_bw()->fig2

# 
# swc_df %>% 
#   na.omit() %>%
#   group_by(si_code) %>% 
#   summarise(biome = unique(biome),
#             "Pearson's correlation" = cor(swc,swc_shallow_mean)
#             ) %>% 
#   mutate(Biome = recode(biome,
#                         `wood` = 'WOOD',
#                         `temp` = "TEMP",
#                         `bor` = "BOR",
#                         `tro` = "TROP",
#                         `dry` = 'DRY'))->p_suma
# 
# p_suma$Biome <- factor(p_suma$Biome,
#                        levels = c("DRY",
#                                   "WOOD",
#                                   "TEMP",
#                                   "BOR",
#                                   "TROP"))
# 
# library("reshape2")
# mm <- melt(p_suma[,3:4])
# ggplot(mm,aes(x=value, fill = Biome))+
#   geom_histogram(position="stack", color = "grey30")+
#   geom_vline(data = p_suma, mapping= aes(xintercept=median(`Pearson's correlation`)),
#              color="blue", linetype="dashed", size=1)+
#   xlab("Pearson's correlation")+
#   theme_bw()+
#   scale_fill_manual(values=custom_palette)-> fig2
# 
# 
# prow <- plot_grid(
#   fig1,#+ theme(legend.position="none"),
#   fig2,
#   align = 'vh',
#   labels = c("A", "B"),
#   hjust = -1,
#   nrow = 1
# )
# prow

ggsave(
  "plots/swc_era5.png",
  fig1,
  # prow,
  width = 13,
  height = 10,
  dpi = 600,
  units = "cm"
)
