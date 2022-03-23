###
###
#' 
#' Script for:
#' A global meta-analysis reveals higher life-history phenotypic variation in urban birds than in their non-urban neighbours
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.09.24.461498
#' 
#' Latest update: 2022/03/20
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Script for among-year meta-model of CVR in laying date, cluth size and number of fledglings
#' 
##
##

##
##### libraries #####
##
pacman::p_load(dplyr, metafor, orchaRd) 

##
##### data #####
##
data <- readRDS("./data/processed_RDS_data_files/metaanalysis_AMONGYEAR_data.RDS")
data$obsID <- 1:nrow(data)
data$scientific_name <- gsub(x = data$scientific_name, 
                             pattern = "_", 
                             replacement = " ")
data$scientific_name_phylo <- data$scientific_name
head(data)

# phylogentic correlation
phylo_cor <- readRDS("./data/processed_RDS_data_files/phylogenetic_correlations_AMONG_model.RDS")


##
## data for model
df_lnCVR_among <- data %>% 
  filter(!is.na(lnCVR)) %>% 
  filter(!is.na(lnCVR.sv))

##
##### Meta-analysis lnCVR among seasons #####
##
lnCVR_among <- rma.mv(yi = lnCVR, 
                      V = lnCVR.sv, 
                      mods = ~ trait - 1, 
                      random =  list(~trait|study_ID,
                                     ~trait|obsID, 
                                     ~1|Pop_ID,
                                     ~1|scientific_name_phylo,
                                     ~1|scientific_name),
                      R = list(scientific_name_phylo = phylo_cor), #phylogenetic relatedness
                      struct=c("DIAG", "DIAG"), 
                      data=df_lnCVR_among, 
                      method="REML")
summary(lnCVR_among) # results presented in Table 1
