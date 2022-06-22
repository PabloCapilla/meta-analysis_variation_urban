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
#' Script for meta-regressions of lnCVR including land cover data as moderators 
#' (difference in urban index and habitat heterogeneity).
#' 
##
##

##
##### libraries and help functions #####
##
pacman::p_load(dplyr, tidyr, extrafont, metafor, ggplot2, orchaRd, data.table) 
loadfonts()
source("./scripts/R_library/functions.R")
source("./scripts/R_library/orchard_plot_PCL.R")
source("./scripts/R_library/orchard_plot_PCL_noApples.R")
source("./scripts/R_library/FUNCTION_lc_selection.R")

##
##### data #####
##
data <- readRDS("./data/processed_RDS_data_files/metaanalysis_full_data_landcover.RDS") # data including landcover variables

# matrix with phylogentic correlations
phylo_cor <- readRDS("./data/processed_RDS_data_files/phylogenetic_correlations_lancover_model.RDS")

##
##
##### Sample size (choose one buffer randomly to summarise sample sizes) #####
##
##

## 'lc_selection' function in 'R_library' folder

df_samplesize <- lc_selection (df = data, 
                               buffer_size = 2000) 
range(df_samplesize$year) # year range
nrow(df_samplesize) # number of observations
length(unique(df_samplesize$study_ID)) # number of studies
length(unique(df_samplesize$species_ID)) # number of species

# differences in urbanisation score at buffer 2000m
df_samplesize %>% 
  group_by(study_ID) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  summarise(mean_urban = mean(urban_urban_value),
            se_urban = sd(urban_urban_value)/sqrt(n()),
            mean_nonurban = mean(nonurban_urban_value),
            se_nonurban = sd(nonurban_urban_value)/sqrt(n()),)


##
##
##### Loop to run models with different buffer sizes #####
##
##
buffer_size <- seq(250, 5000,by = 250) # buffer sizes
df00 <- as.list(NA)

for(m in 1:length(buffer_size)){
  # function to select landcover data from each buffer size
  df_analysis <- lc_selection (df = data, 
                               buffer_size = buffer_size[m]) 
  
  
  
  ## clean and accurate data for model
  df_CVR_exact <- df_analysis %>% 
    select(exact_urban_coor, exact_rural_coor, study_ID, Pop_ID, scientific_name_phylo, 
           scientific_name, trait, lnCVR, lnCVR.sv,
           urban_urban_value, urban_het, nonurban_urban_value, nonurban_het) %>% 
    filter(!is.na(lnCVR)) %>% 
    filter(!is.na(lnCVR.sv)) %>%  
    filter(exact_urban_coor == "YES") %>% 
    filter(exact_rural_coor == "YES") %>% 
    mutate(urban_dif = urban_urban_value - nonurban_urban_value,
           het_dif = urban_het - nonurban_het)
  df_CVR_exact$obsID <- 1:nrow(df_CVR_exact)
  
  # difference in urban index across populations (including one observation per study)
  m_urban <- lm(urban_dif ~ 1, 
                data = df_CVR_exact %>% 
                  group_by(study_ID) %>% 
                  filter(row_number() == 1))
  m_urban_summary <- summary(m_urban)
  
  # meta-model
  model_metaregression <- rma.mv(yi = lnCVR, 
                                 V = lnCVR.sv, 
                                 mods = ~
                                   trait - 1 +
                                   urban_dif +
                                   het_dif,
                                 random =  list(~trait|study_ID,
                                                ~trait|obsID, 
                                                ~1|Pop_ID,
                                                ~1|scientific_name_phylo,
                                                ~1|scientific_name),
                                 R = list(scientific_name_phylo = phylo_cor), 
                                 struct=c("DIAG", "DIAG"), 
                                 data=df_CVR_exact, 
                                 method="REML", 
                                 control=list(optimizer="BFGS", iter.max=1000, rel.tol=1e-8)) # to improve convergence
  
  # saving model output
  file_name <- paste0("LC_model_Buffer", buffer_size[m], ".RDS")
  saveRDS(model_metaregression, file = paste0("./models/land_cover_models/", file_name))
  model_metaregression_est <- estimates.CI(model_metaregression)
  
  # save model results
  df00[[m]] <- model_metaregression_est
  
  df00[[m]]["dif_urban"] <- m_urban_summary$coefficients[1,1]
  df00[[m]]["se_dif_urban"] <- m_urban_summary$coefficients[1,2]
  df00[[m]]["p_dif_urban"] <- m_urban_summary$coefficients[1,4]
  df00[[m]]["buffer"] <- buffer_size[m]
  
  # message
  cat(paste0("Model for buffer ", buffer_size[m], " finished", "\n"))
  
}
df <- rbindlist(df00)

## save results
#saveRDS(object = df, 
 #      file = "./models/land_cover_models/LC_model_summary_lnCVR.RDS")


##
## read results
df <- readRDS("./models/land_cover_models/LC_model_summary_lnCVR.RDS")
head(df)

##
##
##### Plotting results #####
##
##

## diff in urban index
labels_x <- seq(250, 5000, 250)
labels_x[seq(2,20,2)] <- " "

##
## Plot for Figure 4a - difference in urban index between populations
urban_dif_plot <- ggplot(data = df %>% 
                           group_by(buffer) %>% 
                           filter(row_number() == 1), 
                         aes(x = buffer, y = dif_urban, fill = estimate)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  geom_linerange(aes(ymin = dif_urban - se_dif_urban, ymax = dif_urban + se_dif_urban),
                 position = position_dodge(width = 25)) +  
  geom_point(size = 4, shape = 21, 
             position = position_dodge(width = 25),
             color = "black") +
  scale_x_continuous(breaks = seq(250, 5000, 250), 
                     labels = labels_x) +
  scale_fill_manual(values = "#8DA0CB") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title = element_text(size = 14,family = "Arial"),
        panel.grid = element_blank()) +
  labs(x = expression(atop("Spatial scale",
                           "(buffer radius [m])")), 
       y = expression(atop("Difference in urban index between", 
                           "urban and non-urban populations (± SE)")))

# saving plot
ggsave(filename = "./plots/Figure 4a.jpeg", 
       plot = urban_dif_plot, 
       device = "jpeg", 
       width = 80, 
       height = 120, 
       units = "mm")

##
## Plot for Figure 4b - effect of urban index on lnCVR
urban_effect <- ggplot(data = df %>% 
                         filter(estimate == "urban_dif"), 
                       aes(x = buffer, y = mean, fill = estimate)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 position = position_dodge(width = 25)) +  
  geom_point(size = 4, shape = 21, 
             position = position_dodge(width = 25),
             color = "black") +
  scale_y_continuous(limits = c(-0.7, 0.7)) +
  scale_x_continuous(breaks = seq(250, 5000, 250), 
                     labels = labels_x) +
  scale_fill_manual(values = "#8DA0CB") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title = element_text(size = 14,family = "Arial"),
        panel.grid = element_blank()) +
  labs(x = expression(atop("Spatial scale",
                           "(buffer radius [m])")),
       y = expression(atop("Effect of differences in urban index", 
                           "on phenotypic lnCVR (± 95% CI)")))

ggsave(filename = "./plots/Figure 4b.jpeg", 
       plot = urban_effect, 
       device = "jpeg", 
       width = 80, 
       height = 120, 
       units = "mm")

##
## Plot for Figure 4c - effect of difference in heterogeneit on lnCVR
het_effect <- ggplot(data = df %>% 
                       filter(estimate == "het_dif"), 
                     aes(x = buffer, y = mean, fill = estimate)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 position = position_dodge(width = 25)) +  
  geom_point(size = 4, shape = 21, 
             position = position_dodge(width = 25),
             color = "black") +
  scale_y_continuous(limits = c(-0.13, 0.07)) +
  scale_x_continuous(breaks = seq(250, 5000, 250), 
                     labels = labels_x) +
  scale_fill_manual(values = "#8DA0CB") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title = element_text(size = 14,family = "Arial"),
        panel.grid = element_blank()) +
  labs(x = expression(atop("Spatial scale",
                           "(buffer radius [m])")), 
       y = expression(atop("Effect of differences in habitat", 
                           "heterogeneity on lnCVR (± 95% CI)")))

# saving plot
ggsave(filename = "./plots/Figure 4c.jpeg", 
       plot = het_effect, 
       device = "jpeg", 
       width = 80, 
       height = 120, 
       units = "mm")



