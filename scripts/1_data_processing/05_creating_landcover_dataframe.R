###
###
#' 
#' Script for:
#' A global meta-analysis reveals higher phenotypic variation in urban birds than in their non-urban neighbours
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.09.24.461498
#' 
#' Latest update: 2022/06/20
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' This script prepares the dataset for analysis of land cover on lnCVR
#' 
##
##
##


##
##### libraries #####
##
pacman::p_load(tidyverse, extrafont, metafor, orchaRd, RColorBrewer, data.table) 
loadfonts()
source("./scripts/R_library/functions.R")
source("./scripts/R_library/orchard_plot_PCL.R")
source("./scripts/R_library/orchard_plot_PCL_noApples.R")


##
##
##### ES data set #####
##
##
data <- readRDS("./data/processed_RDS_data_files/metaanalysis_full_data.RDS")

##
##
##### LC data #####
##
##
load(file = "./data/processed_RDS_data_files/land_cover_data.RDA")
df_lc <- final_data
head(df_lc)

## reformat dataset for analysis
results_per_study_year <- as.list(NA)


for(i in 1:length(unique(data$study_ID))){
  study_id_loop <- unique(data$study_ID)[i]
  df_loop <- data[data$study_ID == study_id_loop,]
  
  # is there land cover data for study_id_loop?
  if(nrow(df_lc[df_lc$study_ID == study_id_loop,]) == 0){
    next()
  } else {
    
    
    # is that a multiple year study?
    if(any(df_loop$multiple_year_ES == "YES")){
      
      results_per_study_year[[i]] <- left_join(x = data %>% 
                                                 mutate(year = as.character(year)) %>% 
                                                 filter(study_ID == study_id_loop),
                                               y = df_lc %>% 
                                                 filter(study_ID == study_id_loop) %>% 
                                                 group_by(study_ID, buffer) %>% 
                                                 summarise(urban_urbanization  = mean(urban_urbanization, na.rm = T),
                                                           urban_heterogeneity   = mean(urban_heterogeneity , na.rm = T),
                                                           nonurban_urbanization   = mean(nonurban_urbanization, na.rm = T),
                                                           nonurban_heterogeneity  = mean(nonurban_heterogeneity, na.rm = T)) %>% 
                                                 pivot_wider(id_cols = c(study_ID), 
                                                             names_from = buffer, 
                                                             values_from = 3:6),
                                               by = "study_ID") 
      
    } else {
      results_per_study_year[[i]] <- left_join(x = data %>% 
                                                 mutate(year = as.character(year)) %>% 
                                                 filter(study_ID == study_id_loop),
                                               y = df_lc %>% 
                                                 mutate(year = as.character(year)) %>% 
                                                 filter(study_ID == study_id_loop) %>% 
                                                 group_by(study_ID, year, buffer) %>% 
                                                 summarise(urban_urbanization  = mean(urban_urbanization, na.rm = T),
                                                           urban_heterogeneity   = mean(urban_heterogeneity , na.rm = T),
                                                           nonurban_urbanization   = mean(nonurban_urbanization, na.rm = T),
                                                           nonurban_heterogeneity  = mean(nonurban_heterogeneity, na.rm = T)) %>% 
                                                 pivot_wider(id_cols = c(study_ID, year), 
                                                             names_from = buffer, 
                                                             values_from = 4:7),
                                               by = c("study_ID", "year")) 
    }
  }
  
  # control progress
  cat(round(i/nrow(data), digits = 2))
}


# final dataset for landcover analysis
df_final <- results_per_study_year %>% 
  discard(is.null)
df_final <- df_final[-1]
data_lc_final <- rbindlist(df_final)


##
## save dataset
saveRDS(object = data_lc_final, file = "./data/processed_RDS_data_files/metaanalysis_full_data_landcover.RDS")
