###
###
#' 
#' Script for:
#' A global meta-analysis reveals higher variation in breeding phenology in urban birds than in their non-urban neighbours
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.09.24.461498
#' 
#' Latest update: 2022/07/26
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())
##
##
##### Script aim: #####
#' Sensitivity analysis of lnRR and lnCVR results.
#' * lnRR analysis repeated using SMDH
#' * lnCVR analysis repeated using lnVR and arm-based model following Senior et al. 2016 Evolution, Medicine and Public Health
#' 
##
##

##
##### libraries #####
##
pacman::p_load(dplyr, tidyr, extrafont, metafor, ggplot2, orchaRd, DT, RColorBrewer) 
loadfonts()
source("./scripts/R_library/functions.R")
source("./scripts/R_library/orchard_plot_PCL.R")
source("./scripts/R_library/orchard_plot_PCL_noApples.R")

##
##
##### Data #####
##
##
# ES data
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
##
##### Sensitivity analysis of lnRR #####
##
##

## data for lnRR
df_SMDH <- data %>%
  filter(!is.na(SMDH)) %>% 
  filter(!is.na(SMDH.sv))

# SMDH model (same model that used for lnRR - MODEL 2)
model_mean_SMDH <- rma.mv(yi = SMDH, 
                 V = SMDH.sv, 
                 mods = ~ trait - 1, 
                 random =  list(~trait|study_ID,
                                ~trait|obsID, 
                                ~1|Pop_ID,
                                ~1|scientific_name_phylo,
                                ~1|scientific_name),
                 R = list(scientific_name_phylo = phylo_cor), #phylogenetic relatedness
                 struct=c("UN", "DIAG"), 
                 data=df_SMDH,
                 method="REML")
summary(model_mean_SMDH)

orchard_plot_PCL(object = model_mean_SMDH, 
                 mod = " ", 
                 est_point_size = 5,
                 alpha = 0.5,
                 cb = FALSE,
                 xlab = "SMDH",
                 ylab = "Intercept",
                 transfm = "none") +
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 12),
        axis.text.y = element_text("Arial", size = 12, angle = 45),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = rev(brewer.pal(n = 3, "Set2"))) +
  scale_color_manual(values = rev(brewer.pal(n = 3, "Set2"))) +
  scale_y_discrete(labels = c("Laying date", "Clutch size", "# fledlgings"))


##
##
##### Sensitivity analysis of lnCVR #####
##
##

## data for lnVR
df_lnVR <- data %>%
  filter(!is.na(lnVR)) %>% 
  filter(!is.na(lnVR.sv))

##
## Analysis of lnVR
model_var_lnVR <- rma.mv(yi = lnVR, 
                         V = lnVR.sv, 
                         mods = ~ trait - 1, 
                         random =  list(~trait|study_ID,
                                        ~trait|obsID, 
                                        ~1|Pop_ID,
                                        ~1|scientific_name_phylo,
                                        ~1|scientific_name),
                         R = list(scientific_name_phylo = phylo_cor), #phylogenetic relatedness
                         struct=c("DIAG", "DIAG"), 
                         data=df_lnVR, 
                         method="REML")
summary(model_var_lnVR)

## 
## plot for lnVR
orchard_plot_PCL(object = model_var_lnVR, 
                 mod = " ", 
                 est_point_size = 5,
                 alpha = 0.5,
                 cb = FALSE,
                 xlab = "log Total Variation Ratio (lnVR)",
                 ylab = "Intercept",
                 transfm = "none") +
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 12),
        axis.text.y = element_text("Arial", size = 12, angle = 45),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = rev(brewer.pal(n = 3, "Set2"))) +
  scale_color_manual(values = rev(brewer.pal(n = 3, "Set2"))) +
  scale_y_discrete(labels = c("Laying date", "Clutch size", "# fledlgings"))



##
## Arm-based model ln lnSD

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
## Arm-based model for lnSD in laying date
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

##
## Arm-based model for lnSD in clutch size
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

##
## Arm-based model for lnSD in fledlging number
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




