###
###
#' 
#' Script for:
#' A global meta-analysis reveals higher life-history phenotypic variation in urban birds than in their non-urban neighbours
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.09.24.461498
#' 
#' Latest update: 2022/03/18
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Meta-analysis of lnRR (presented in section i of the manuscript)
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
## data for lnRR
df_lnRR <- data %>%
  filter(!is.na(lnRR)) %>% 
  filter(!is.na(lnRR.sv))

# manually inspect very large lnRR 
df_lnRR %>% 
  filter(lnRR < -2) # manually inspected and seemingly fine


##
##### Meta-analysis lnRR #####
##
model_lnRR <- rma.mv(yi = lnRR, 
                      V = lnRR.sv, 
                      mods = ~ 1, 
                      random = list(~1|study_ID,
                                    ~1|Pop_ID,
                                    ~1|obsID, 
                                    ~1|scientific_name_phylo,
                                    ~1|scientific_name),
                      R = list(scientific_name_phylo = phylo_cor),
                      data=df_lnRR, 
                      method="REML")

summary(model_lnRR) # model summary
round(i2_ml(model_lnRR), digits = 4) # heterogeneity random effects

##
##### Meta-model lnRR without imputed values #####
##

#
# removing imputed SD values from dataset produces virtually same results
df_lnRR_no_input <- df_lnRR %>% 
  filter(sd_imputation_urban == "NONE") %>% 
  filter(sd_imputation_rural == "NONE")

model_lnRR_no_input <- rma.mv(yi = lnRR, 
                         V = lnRR.sv, 
                         mods = ~ 1, 
                         random = list(~1|study_ID,
                                       ~1|Pop_ID,
                                       ~1|obsID, 
                                       ~1|scientific_name_phylo,
                                       ~1|scientific_name),
                         R = list(scientific_name_phylo = phylo_cor),
                         data=df_lnRR_no_input, 
                         method="REML") # model produces warning as some species in the matrix don't appear in trimmed dataset
summary(model_lnRR_no_input) # model summary
round(i2_ml(model_lnRR_no_input), digits = 4) # heterogeneity random effects

##
##
##### Plotting results #####
##
##

## plot lnRR overall
lnRR_plot <- orchard_plot_PCL(object = model_lnRR, 
                              mod = "Intercept", 
                              est_point_size = 5,
                              alpha = 0.5,
                              pred_int = T,
                              cb = FALSE,
                              xlab = "log Response Ratio", 
                              ylab = "Intercept",
                              transfm = "none") +
  scale_fill_manual(values = "#bdbdbd") +
  scale_color_manual(values = "#bdbdbd") 

lnRR_plot <- lnRR_plot + 
  theme(axis.title = element_text("Arial", size = 10),
        axis.text.x = element_text("Arial", size = 10),
        axis.text.y = element_blank(),
        legend.position = "none") 

ggsave(filename = "./plots/Figure S4.jpeg", 
       plot = lnRR_plot, 
       device = "jpeg", 
       height = 90, 
       width = 150, 
       units = "mm")

# zoom in of -0.5 < x < 0.5
lnRR_plot_close <- lnRR_plot + 
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none") +
  scale_x_continuous(limits = c(-0.50, 0.50)) 
  

ggsave(filename = "./plots/Figure S4 snippet.jpeg", 
       plot = lnRR_plot_close, 
       device = "jpeg", 
       height = 50, 
       width = 100, 
       units = "mm")

