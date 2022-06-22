###
###
#' 
#' Script for:
#' A global meta-analysis reveals higher phenological variation in urban birds than in their non-urban neighbours
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.09.24.461498
#' 
#' Latest update: 2022/06/22
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Meta-regressions to assess publication bias in lnRR (presented in section iv of the manuscript)
#' 
##
##

##
##### libraries #####
##
pacman::p_load(dplyr, tidyr, metafor, orchaRd) 

##
##### data #####
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
## data for model
df_lnRR_pb <- data %>% 
  filter(!is.na(lnRR)) %>% 
  filter(!is.na(lnRR.sv)) %>% 
  mutate(s_weights = 1/lnRR.sv,
         sqrt_inv_eff_ss = sqrt((1/urban_n) + (1/rural_n)),
         inv_eff_ss = (1/urban_n) + (1/rural_n),
         publication_year_c = publication_year - mean(publication_year))

##
##
##### Publication bias #####
##
##

## small-study effects
model_pub_bias1 <- rma.mv(yi = lnRR, 
                         V = lnRR.sv, 
                         mods = ~ 1 +
                           sqrt_inv_eff_ss, 
                         random = list(~1|study_ID,
                                       ~1|Pop_ID,
                                       ~1|obsID, 
                                       ~1|scientific_name_phylo,
                                       ~1|scientific_name),
                         R = list(scientific_name_phylo = phylo_cor),
                         data=df_lnRR_pb, 
                         method="REML")
summary(model_pub_bias1)
r2_ml(model_pub_bias1)


## time lag effects
model_pub_bias2 <- rma.mv(yi = lnRR, 
                          V = lnRR.sv, 
                          mods = ~ 1 +
                            publication_year_c, 
                          random = list(~1|study_ID,
                                        ~1|Pop_ID,
                                        ~1|obsID, 
                                        ~1|scientific_name_phylo,
                                        ~1|scientific_name),
                          R = list(scientific_name_phylo = phylo_cor),
                          data=df_lnRR_pb, 
                          method="REML")
summary(model_pub_bias2)
r2_ml(model_pub_bias2)

