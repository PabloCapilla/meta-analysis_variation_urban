###
###
#' 
#' Script for:
#' A global meta-analysis reveals higher life-history phenotypic variation in urban birds than in their non-urban neighbours
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.09.24.461498
#' 
#' Latest update: 2022/06/16
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' This script creates a dataframe to investigate lnCVR within breeding seasons
#' 
##
##
##

##
##### libraries #####
##
pacman::p_load(dplyr, tidyr) 
loadfonts()

##
##### data #####
##

# Full data set
data <- readRDS("./data/processed_RDS_data_files/metaanalysis_full_data.RDS")

# filter data from individual years
single_year <- data %>% filter(multiple_year_ES == "NO")

# save table
saveRDS(object = single_year, file = "./data/processed_RDS_data_files/metaanalysis_WITHINYEAR_data.RDS")
