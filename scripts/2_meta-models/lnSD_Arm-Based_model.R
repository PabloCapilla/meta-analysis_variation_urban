###
###
#' 
#' Script for:
#' A global meta-analysis reveals higher life-history phenotypic variation in urban birds than in their non-urban neighbours
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.09.24.461498
#' 
#' Latest update: 2022/0608
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Meta-analysis of lnSD using arm-based model 
#' (see Senior et al. 2016 and Eq. 18 in Nakagawa et al. 2015). 
#' This additional model is used as a sensitivity analysis to confirm our results analysing lnCVR.
#' 
##
##

##
##### libraries #####
##
pacman::p_load(dplyr, tidyr, extrafont, metafor, ggplot2, orchaRd) 
loadfonts()
source("./scripts/R_library/orchard_plot_PCL.R")

##
##### data #####
##
data <- readRDS("./data/processed_RDS_data_files/metaanalysis_full_data.RDS")
data$obsID <- 1:nrow(data)
data$scientific_name <- gsub(x = data$scientific_name, 
                             pattern = "_", 
                             replacement = " ")
data$scientific_name_phylo <- data$scientific_name
head(data)

# matrix with phylogentic correlations
phylo_cor <- readRDS("./data/processed_RDS_data_files/phylogenetic_correlations.RDS")

##
## data for model
## (our original data frame needs to be re-formated to run this analysis)

# urban dataset
dataSD_urban <- data %>% 
  select(study_ID, Pop_ID, scientific_name_phylo, scientific_name, trait, year,
         urban_mean, urban_sd, urban_n) %>% 
  mutate(habitat = "urban") %>%
  rename(mean = urban_mean, sd = urban_sd, n = urban_n)

# rural dataset
dataSD_rural <- data %>% 
  select(study_ID, Pop_ID, scientific_name_phylo, scientific_name, trait, year,
         rural_mean, rural_sd, rural_n) %>% 
  mutate(habitat = "rural") %>%
  rename(mean = rural_mean, sd = rural_sd, n = rural_n)

# re-formated dataset
dataSD <- rbind(dataSD_urban, dataSD_rural)

##
## calculation of lnSD and sampling variance in lnSD (following Eq8 in Nakagawa et al 2015)
dataSD <- dataSD %>% 
  mutate(lnSD = log(sd) + 1 / (2*(n-1)),
         sampling_variance_SD = 1 / (2*(n-1)), 
         lnMEAN = log(mean),
         obsID = 1:nrow(dataSD))


##
##### Arm-based model for lnSD in laying date #####
##
model_layingdate_lnSD <- rma.mv(yi = lnSD, 
                                V = sampling_variance_SD, 
                                mods = ~ habitat + lnMEAN + 1, 
                                random = list(~habitat|study_ID,
                                              ~1|Pop_ID,
                                              ~1|obsID, 
                                              ~1|scientific_name_phylo,
                                              ~1|scientific_name),
                                R = list(scientific_name_phylo = phylo_cor), 
                                data = dataSD %>% 
                                  filter(trait == "Laying date"),
                                method = "REML")


summary(model_layingdate_lnSD) 
round(i2_ml(model_layingdate_lnSD), digits = 4)
round(r2_ml(model_layingdate_lnSD), digits = 4)

model_clutchsize_lnSD <- rma.mv(yi = lnSD, 
                                V = sampling_variance_SD, 
                                mods = ~ habitat + lnMEAN + 1, 
                                random = list(~habitat|study_ID,
                                              ~1|Pop_ID,
                                              ~1|obsID, 
                                              ~1|scientific_name_phylo,
                                              ~1|scientific_name),
                                R = list(scientific_name_phylo = phylo_cor), 
                                data = dataSD %>% 
                                  filter(trait == "Clutch size"),
                                method = "REML")
summary(model_clutchsize_lnSD) 
round(i2_ml(model_lnCVR), digits = 4)


model_fledglings_lnSD <- rma.mv(yi = lnSD, 
                                V = sampling_variance_SD, 
                                mods = ~ habitat + lnMEAN + 1, 
                                random = list(~habitat|study_ID,
                                              ~1|Pop_ID,
                                              ~1|obsID, 
                                              ~1|scientific_name_phylo,
                                              ~1|scientific_name),
                                R = list(scientific_name_phylo = phylo_cor), 
                                data = dataSD %>% 
                                  filter(trait == "# Fledglings"),
                                method = "REML")
summary(model_fledglings_lnSD) 
round(i2_ml(model_fledglings_lnSD), digits = 4)


