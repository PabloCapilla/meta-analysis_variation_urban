###
###
#' 
#' Script for:
#' A global meta-analysis reveals higher life-history phenotypic variation in urban birds than in their non-urban neighbours
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.09.24.461498
#' 
#' Latest update: 2022/03/19
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' This script creates a dataframe to investigate lnCVR among breeding seasons
#' 
##
##
##

##
##### libraries & functions #####
##
pacman::p_load(dplyr, tidyr, metafor, ggplot2, orchaRd, stringr) 
source("./scripts/R_library/FUNCTION_combined_groups.R")

##
##### data #####
##

# Full data set
data <- readRDS("./data/processed_RDS_data_files/metaanalysis_full_data.RDS")

# store multi-year and single year data
multi_year <- data %>% filter(multiple_year_ES == "YES")
single_year <- data %>% filter(multiple_year_ES == "NO")

##
##### summary of the number of years per study (in those with single-year estimates) #####
##
sum_studies <- single_year %>% 
  group_by(study_ID, trait, species_ID) %>% 
  summarise(n_per_study = n()) %>% 
  arrange(desc(n_per_study)) %>% 
  filter(n_per_study > 1) 
table(sum_studies$trait)


##
#### Function to calculate mean and sd among years from single year data #####
##
study_ids <- unique(sum_studies$study_ID)
store_list <- list(NA) # list to store new data

# loop to calculate means and sd
for(s in 1:length(study_ids)){
  # study for loop s
  study_loop <- study_ids[s]
  
  # data from study s
  df_loop_study <- single_year[single_year$study_ID == study_loop, ]
  df_loop_PopID <- unique(df_loop_study$Pop_ID)
  trait_ids <- unique(df_loop_study$trait)
  
  for(t in 1:length(trait_ids)){
    
    # trait for loop t
    trait_loop <- trait_ids[t]
    
    # data from study s
    df_loop <- single_year[single_year$study_ID == study_loop & single_year$trait == trait_loop, ]
    
    # are there data in 'df_loop'?
    if(length(df_loop) == 0){
      next()
    } else {
      # are there data for more than one species?
      for(b in 1:length(unique(df_loop$species_ID))){
        df_calc <- df_loop[df_loop$species_ID == unique(df_loop$species_ID)[b],]
        
        # years
        year_vector <- paste0(as.numeric(df_calc$year), collapse = "_")
        
        # urban combined estimates
        urban_combined <- combine_groups(means = df_calc$urban_mean, 
                                         sds = df_calc$urban_sd, 
                                         ns = df_calc$urban_n)
        urban_combined_n <- sum(df_calc$urban_n)
        
        # rural combined estimates
        rural_combined <- combine_groups(means = df_calc$rural_mean, 
                                         sds = df_calc$rural_sd, 
                                         ns = df_calc$rural_n)
        rural_combined_n <- sum(df_calc$rural_n)
        
        # data to save
        length_list <- length(store_list)
        store_list[[length_list + 1]] <- data.frame(study_ID = study_loop,
                                                    species_ID = unique(df_loop$species_ID)[b],
                                                    scientific_name = unique(df_calc$scientific_name),
                                                    Pop_ID = df_loop_PopID,
                                                    trait = trait_loop,
                                                    year = year_vector,
                                                    n_year = length(df_calc$year),
                                                    multiple_year_ES = "YES",
                                                    urban_mean = as.numeric(urban_combined[1]),
                                                    urban_sd = as.numeric(urban_combined[2]),
                                                    urban_n = urban_combined_n,
                                                    rural_mean = as.numeric(rural_combined[1]),
                                                    rural_sd = as.numeric(rural_combined[2]),
                                                    rural_n = rural_combined_n,
                                                    first_year = unique(df_calc$first_year),
                                                    publication_year = unique(df_calc$publication_year),
                                                    urban_lat = unique(df_calc$urban_lat),
                                                    urban_lon = unique(df_calc$urban_lon),
                                                    exact_urban_coor = unique(df_calc$exact_urban_coor),
                                                    rural_lat = unique(df_calc$rural_lat),
                                                    rural_lon = unique(df_calc$rural_lon),
                                                    exact_rural_coor = unique(df_calc$exact_rural_coor))
      }
    }
  }
}
among_df <- do.call(what = "rbind", args = store_list[!is.na(store_list)])


##
##### Calculating meta-analysis y variables #####
##

## lnCVR
among_df <- as.data.frame(escalc(measure="CVR", 
                             n1i=urban_n, n2i=rural_n, 
                             m1i=urban_mean, m2i=rural_mean, 
                             sd1i=urban_sd, sd2i=rural_sd,
                             data=among_df,
                             var.names=c("lnCVR","lnCVR.sv"), 
                             add.measure=F,
                             append=TRUE))

## lnVR
among_df <- as.data.frame(escalc(measure="VR", 
                             n1i=urban_n, n2i=rural_n, 
                             m1i=urban_mean, m2i=rural_mean, 
                             sd1i=urban_sd, sd2i=rural_sd,
                             data=among_df,
                             var.names=c("lnVR","lnVR.sv"), 
                             add.measure=F,
                             append=TRUE))

## lnRR
among_df <- as.data.frame(escalc(measure="ROM", 
                             n1i=urban_n, n2i=rural_n, 
                             m1i=urban_mean, m2i=rural_mean, 
                             sd1i=urban_sd, sd2i=rural_sd,
                             data=among_df,
                             var.names=c("lnRR","lnRR.sv"), 
                             add.measure=F,
                             append=TRUE))
head(among_df)

# flag to signal original data was in 'single-year' format
among_df$method_combine <- "calculated"


# combining new multi-year data and original data with multi year estimations
multi_year$n_year <- do.call(what = "c",
                             args = lapply(X = str_split(string = multi_year$year, pattern = "_"), 
                                           FUN = length))

multi_year <- multi_year %>% 
  select(study_ID, species_ID, scientific_name, Pop_ID, trait, year, n_year, multiple_year_ES,
         urban_mean, urban_sd, urban_n, rural_mean, rural_sd, rural_n,
         first_year, publication_year,
         urban_lat, urban_lon, exact_urban_coor, rural_lat, rural_lon, exact_rural_coor,
         lnCVR, lnCVR.sv, lnVR, lnVR.sv, lnRR, lnRR.sv) %>% 
  mutate(method_combine = "original")

data_combined <- rbind(among_df,multi_year)
                          

## quick visualisation of new data set
ggplot(data = data_combined, aes(x = trait, y = lnCVR, fill = method_combine)) +
  geom_boxplot() +
  theme_bw()

ggplot(data = data_combined, aes(x = trait, y = lnRR, lnCVR, fill = method_combine)) +
  geom_boxplot() +
  theme_bw()

ggplot(data = data_combined, aes(x = trait, y = lnVR, lnCVR, fill = method_combine)) +
  geom_boxplot() +
  theme_bw()



##
## save full table
saveRDS(object = data_combined, file = "./data/processed_RDS_data_files/metaanalysis_AMONGYEAR_data.RDS")

