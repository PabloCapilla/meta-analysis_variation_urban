###
###
#' 
#' Script for:
#' A global meta-analysis reveals higher phenological variation in urban birds than in their non-urban neighbours
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.09.24.461498
#' 
#' Latest update: 2022/06/21
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' This script creates a matrix of phylogenetic correlation to include in meta-analytical models of lnCVR including landcover moderators
#' 
##
##
##


##
##### libraries #####
##
pacman::p_load(dplyr, tidyr, metafor, rotl, 
               ape, curl, treeio,
               phytools)

##
##### data #####
##
data <- readRDS("./data/processed_RDS_data_files/landcover_dataframe.RDS")
data <- data %>% 
  filter(!is.na(lnCVR)) %>% 
  filter(!is.na(lnCVR.sv)) %>%  
  filter(exact_urban_coor == "YES") %>% 
  filter(exact_rural_coor == "YES") 
head(data)

##
##### recovering phylogenetic relationship #####
##
data$scientific_name <- gsub(x = data$scientific_name, 
                             pattern = "_", 
                             replacement = " ")
taxa.corrected <- tnrs_match_names(names = unique(data$scientific_name))

# check approximate matches: OK
taxa.corrected[taxa.corrected$approximate_match==TRUE,]

# check synonyms matches: OK
taxa.corrected[taxa.corrected$is_synonym==TRUE,]

# check number of matches: OK
taxa.corrected[taxa.corrected$number_matches>1,]

# retrieving phylogenetic relationships among taxa in the form of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa.corrected[["ott_id"]], label_format = "name")

##
##
##### Creating phylogenetic relationship matrix for meta-analyses #####
##
##

# check for the existence of polytomies
is.binary.tree(tree) # there are no polytomies


# checking that the tree includes every species in data table
tree$tip.label <- gsub("_"," ", tree$tip.label)
intersect(as.character(tree$tip.label), as.character(data$scientific_name))
setdiff(as.character(data$scientific_name), as.character(tree$tip.label)) #listed in our database but not in the tree
setdiff(as.character(tree$tip.label),as.character(data$scientific_name)) # listed in the tree but not in our database

# compute branch lengths of tree
phylo_branch <- compute.brlen(tree, method = "Grafen", power = 1)

# check tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

# matrix to be included in the models
phylo_cor <- vcv(phylo_branch, cor = T)

# 
# save matrix for analyses
saveRDS(phylo_cor, file = "./data/processed_RDS_data_files/phylogenetic_correlations_landcover_models.RDS")



