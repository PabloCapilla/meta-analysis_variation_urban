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
#' Meta-analysis of lnCVR (presented in section ii of the manuscript)
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
df_lnCVR <- data %>% 
  filter(!is.na(lnCVR)) %>% 
  filter(!is.na(lnCVR.sv))

##
##### Meta-analysis lnCVR #####
##
model_lnCVR <- rma.mv(yi = lnCVR, 
                      V = lnCVR.sv, 
                      mods = ~1, 
                      random = list(~1|study_ID,
                                    ~1|Pop_ID,
                                    ~1|obsID, 
                                    ~1|scientific_name_phylo,
                                    ~1|scientific_name),
                      R = list(scientific_name_phylo = phylo_cor), 
                      data = df_lnCVR,
                      method = "REML")
summary(model_lnCVR) 
round(i2_ml(model_lnCVR), digits = 4)


##
##
##### Plotting results #####
##
##

## plot lnCVR overall
lnCVR_plot <- orchard_plot_PCL(object = model_lnCVR, 
                              mod = "Intercept", 
                              est_point_size = 5,
                              alpha = 0.5,
                              cb = FALSE,
                              xlab = "log Coefficient of Variation Ratio", 
                              ylab = "Intercept",
                              transfm = "none") +
  scale_fill_manual(values = "#bdbdbd") +
  scale_color_manual(values = "#bdbdbd") 

lnCVR_plot <- lnCVR_plot + 
  theme(axis.title = element_text("Arial", size = 10),
        axis.text = element_text("Arial", size = 10),
        axis.text.y = element_blank(),
        legend.position = "none") 


ggsave(filename = "./plots/Figure S5.jpeg", 
       plot = lnCVR_plot, 
       device = "jpeg", 
       height = 90, 
       width = 150, 
       units = "mm")


