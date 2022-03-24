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
#' Meta-analysis of lnCVR (presented in section ii of the manuscript)
#' Model selection to define the best structure for correlated random effects
#' 
##
##

##
##### libraries and help functions #####
##
pacman::p_load(dplyr, tidyr, extrafont, metafor, ggplot2, orchaRd, DT) 
loadfonts()
source("./scripts/R_library/functions.R")
source("./scripts/R_library/orchard_plot_PCL.R")
source("./scripts/R_library/orchard_plot_PCL_noApples.R")

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
## data for models
df_lnCVR <- data %>% 
  filter(!is.na(lnCVR)) %>% 
  filter(!is.na(lnCVR.sv))

##
##
##### table to store model results #####
##
##
n_models <- 5
df00 <- data.frame(model_ID = rep(NA,n_models),
                   model_ID_n = rep(NA, n_models),
                   study_id_str = rep(NA, n_models),
                   loglik = rep(NA,n_models),
                   df = rep(NA,n_models),
                   AIC = rep(NA,n_models),
                   coef_LD = rep(NA,n_models),
                   coef_LD_low = rep(NA,n_models),
                   coef_LD_high = rep(NA,n_models),
                   coef_CS = rep(NA,n_models),
                   coef_CS_low = rep(NA,n_models),
                   coef_CS_high = rep(NA,n_models),
                   coef_FN = rep(NA,n_models),
                   coef_FN_low = rep(NA,n_models),
                   coef_FN_high = rep(NA,n_models))

##
##
##### Model 0 - Univariate model not accounting for heterogeneous variances #####
##
##
model0 <- rma.mv(yi = lnCVR, 
                 V = lnCVR.sv, 
                 mods = ~ trait - 1, 
                 random =  list(~1|study_ID,
                                ~trait|obsID, 
                                ~1|Pop_ID,
                                ~1|scientific_name_phylo,
                                ~1|scientific_name),
                 R = list(scientific_name_phylo = phylo_cor),
                 struct=c("DIAG"),  
                 data=df_lnCVR, 
                 method="ML")
#saveRDS(object = model0, "./models/Table_S4/lnCVR_model0.RDS")
#model0 <- readRDS("./models/Table_S4/lnCVR_model0.RDS")

summary(model0)
model0_est <- estimates.CI(model0)

# save results model 0
df00$model_ID[1] <- "NONE"
df00$model_ID_n[1] <- 0
df00$study_id_str[1] <- "NONE"
df00$loglik[1] <- as.numeric(logLik(model0))
df00$df[1] <- attr(logLik(model0), "df")
df00$AIC[1] <- AIC(model0)
df00$coef_LD[1] <- model0_est$mean[model0_est$estimate == "traitLaying date"]
df00$coef_LD_low[1] <- model0_est$lower[model0_est$estimate == "traitLaying date"]
df00$coef_LD_high[1] <- model0_est$upper[model0_est$estimate == "traitLaying date"]
df00$coef_CS[1] <- model0_est$mean[model0_est$estimate == "traitClutch size"]
df00$coef_CS_low[1] <- model0_est$lower[model0_est$estimate == "traitClutch size"]
df00$coef_CS_high[1] <- model0_est$upper[model0_est$estimate == "traitClutch size"]
df00$coef_FN[1] <- model0_est$mean[model0_est$estimate == "trait# Fledglings"]
df00$coef_FN_low[1] <- model0_est$lower[model0_est$estimate == "trait# Fledglings"]
df00$coef_FN_high[1] <- model0_est$upper[model0_est$estimate == "trait# Fledglings"]

##
##
##### Model 1 - Initial trivariate model DIAG #####
##
##
model1 <- rma.mv(yi = lnCVR, 
                 V = lnCVR.sv, 
                 mods = ~ trait - 1, 
                 random =  list(~trait|study_ID,
                                ~trait|obsID, 
                                ~1|Pop_ID,
                                ~1|scientific_name_phylo,
                                ~1|scientific_name),
                 R = list(scientific_name_phylo = phylo_cor), #phylogenetic relatedness
                 struct=c("DIAG", "DIAG"), 
                 data=df_lnCVR, 
                 method="ML")
#saveRDS(object = model1, "./models/Table_S4/lnCVR_model1.RDS")
#model1 <- readRDS("./models/Table_S4/lnCVR_model1.RDS")

##
## top model, as shown in Table S4 (see below)
summary(model1)
r2_ml(model1)
model1_est <- estimates.CI(model1)

# save results model 1
df00$model_ID[2] <- "DIAG"
df00$model_ID_n[2] <- 1
df00$study_id_str[2] <- "DIAG"
df00$loglik[2] <- as.numeric(logLik(model1))
df00$df[2] <- attr(logLik(model1), "df")
df00$AIC[2] <- AIC(model1)
df00$coef_LD[2] <- model1_est$mean[model1_est$estimate == "traitLaying date"]
df00$coef_LD_low[2] <- model1_est$lower[model1_est$estimate == "traitLaying date"]
df00$coef_LD_high[2] <- model1_est$upper[model1_est$estimate == "traitLaying date"]
df00$coef_CS[2] <- model1_est$mean[model1_est$estimate == "traitClutch size"]
df00$coef_CS_low[2] <- model1_est$lower[model1_est$estimate == "traitClutch size"]
df00$coef_CS_high[2] <- model1_est$upper[model1_est$estimate == "traitClutch size"]
df00$coef_FN[2] <- model1_est$mean[model1_est$estimate == "trait# Fledglings"]
df00$coef_FN_low[2] <- model1_est$lower[model1_est$estimate == "trait# Fledglings"]
df00$coef_FN_high[2] <- model1_est$upper[model1_est$estimate == "trait# Fledglings"]

##
##
##### Model 2 - Trivariate model CS #####
##
##
model2 <- rma.mv(yi = lnCVR, 
                 V = lnCVR.sv, 
                 mods = ~ trait - 1, 
                 random =  list(~trait|study_ID,
                                ~trait|obsID, 
                                ~1|Pop_ID,
                                ~1|scientific_name_phylo,
                                ~1|scientific_name),
                 R = list(scientific_name_phylo = phylo_cor), #phylogenetic relatedness
                 struct=c("CS", "DIAG"), 
                 data=df_lnCVR, 
                 method="ML")
#saveRDS(object = model2, "./models/Table_S4/lnCVR_model2.RDS")
#model2 <- readRDS("./models/Table_S4/lnCVR_model2.RDS")

summary(model2)
model2_est <- estimates.CI(model2)

# save results model 2
m_row <- 3
df00$model_ID[m_row] <- "CS"
df00$model_ID_n[m_row] <- m_row-1
df00$study_id_str[m_row] <- "CS"
df00$loglik[m_row] <- as.numeric(logLik(model2))
df00$df[m_row] <- attr(logLik(model2), "df")
df00$AIC[m_row] <- AIC(model2)
df00$coef_LD[m_row] <- model2_est$mean[model2_est$estimate == "traitLaying date"]
df00$coef_LD_low[m_row] <- model2_est$lower[model2_est$estimate == "traitLaying date"]
df00$coef_LD_high[m_row] <- model2_est$upper[model2_est$estimate == "traitLaying date"]
df00$coef_CS[m_row] <- model2_est$mean[model2_est$estimate == "traitClutch size"]
df00$coef_CS_low[m_row] <- model2_est$lower[model2_est$estimate == "traitClutch size"]
df00$coef_CS_high[m_row] <- model2_est$upper[model2_est$estimate == "traitClutch size"]
df00$coef_FN[m_row] <- model2_est$mean[model2_est$estimate == "trait# Fledglings"]
df00$coef_FN_low[m_row] <- model2_est$lower[model2_est$estimate == "trait# Fledglings"]
df00$coef_FN_high[m_row] <- model2_est$upper[model2_est$estimate == "trait# Fledglings"]

##
##
##### Model 3 - Trivariate model HCS #####
##
##
model3 <- rma.mv(yi = lnCVR, 
                 V = lnCVR.sv, 
                 mods = ~ trait - 1, 
                 random =  list(~trait|study_ID,
                                ~trait|obsID, 
                                ~1|Pop_ID,
                                ~1|scientific_name_phylo,
                                ~1|scientific_name),
                 R = list(scientific_name_phylo = phylo_cor), #phylogenetic relatedness
                 struct=c("HCS", "DIAG"), 
                 data=df_lnCVR, 
                 method="ML")
#saveRDS(object = model3, "./models/Table_S4/lnCVR_model3.RDS")
#model3 <- readRDS("./models/Table_S4/lnCVR_model3.RDS")

summary(model3)
model3_est <- estimates.CI(model3)

# save results model 3
m_row <- 4
df00$model_ID[m_row] <- "HCS"
df00$model_ID_n[m_row] <- m_row-1
df00$study_id_str[m_row] <- "HCS"
df00$loglik[m_row] <- as.numeric(logLik(model3))
df00$df[m_row] <- attr(logLik(model3), "df")
df00$AIC[m_row] <- AIC(model3)
df00$coef_LD[m_row] <- model3_est$mean[model3_est$estimate == "traitLaying date"]
df00$coef_LD_low[m_row] <- model3_est$lower[model3_est$estimate == "traitLaying date"]
df00$coef_LD_high[m_row] <- model3_est$upper[model3_est$estimate == "traitLaying date"]
df00$coef_CS[m_row] <- model3_est$mean[model3_est$estimate == "traitClutch size"]
df00$coef_CS_low[m_row] <- model3_est$lower[model3_est$estimate == "traitClutch size"]
df00$coef_CS_high[m_row] <- model3_est$upper[model3_est$estimate == "traitClutch size"]
df00$coef_FN[m_row] <- model3_est$mean[model3_est$estimate == "trait# Fledglings"]
df00$coef_FN_low[m_row] <- model3_est$lower[model3_est$estimate == "trait# Fledglings"]
df00$coef_FN_high[m_row] <- model3_est$upper[model3_est$estimate == "trait# Fledglings"]

##
##
##### Model 4 - Trivariate model UN #####
##
##
model4 <- rma.mv(yi = lnCVR, 
                 V = lnCVR.sv, 
                 mods = ~ trait - 1, 
                 random =  list(~trait|study_ID,
                                ~trait|obsID, 
                                ~1|Pop_ID,
                                ~1|scientific_name_phylo,
                                ~1|scientific_name),
                 R = list(scientific_name_phylo = phylo_cor), #phylogenetic relatedness
                 struct=c("UN", "DIAG"), 
                 data=df_lnCVR, 
                 method="ML")
#saveRDS(object = model4, "./models/Table_S4/lnCVR_model4.RDS")
#model4 <- readRDS("./models/Table_S4/lnCVR_model4.RDS")

summary(model4)
model4_est <- estimates.CI(model4)

# save results model 4
m_row <- 5
df00$model_ID[m_row] <- "UN"
df00$model_ID_n[m_row] <- m_row-1
df00$study_id_str[m_row] <- "UN"
df00$loglik[m_row] <- as.numeric(logLik(model4))
df00$df[m_row] <- attr(logLik(model4), "df")
df00$AIC[m_row] <- AIC(model4)
df00$coef_LD[m_row] <- model4_est$mean[model4_est$estimate == "traitLaying date"]
df00$coef_LD_low[m_row] <- model4_est$lower[model4_est$estimate == "traitLaying date"]
df00$coef_LD_high[m_row] <- model4_est$upper[model4_est$estimate == "traitLaying date"]
df00$coef_CS[m_row] <- model4_est$mean[model4_est$estimate == "traitClutch size"]
df00$coef_CS_low[m_row] <- model4_est$lower[model4_est$estimate == "traitClutch size"]
df00$coef_CS_high[m_row] <- model4_est$upper[model4_est$estimate == "traitClutch size"]
df00$coef_FN[m_row] <- model4_est$mean[model4_est$estimate == "trait# Fledglings"]
df00$coef_FN_low[m_row] <- model4_est$lower[model4_est$estimate == "trait# Fledglings"]
df00$coef_FN_high[m_row] <- model4_est$upper[model4_est$estimate == "trait# Fledglings"]

##
##
##### Results #####
##
##
round_to <- 2
min_AIC <- min(df00$AIC)
df00$delta <- round(df00$AIC - min_AIC, digits = round_to)
df00$AIC <- round(df00$AIC, digits = round_to)
df00$loglik <- round(df00$loglik, digits = round_to)

## 
## TABLE S4
##
datatable(
  df00 %>%
    select(`Variance-covariance structure` = model_ID,
           k = df,
           `Log-likelihood` = loglik,
           AIC = AIC,
           dAIC = delta) %>% 
    arrange(dAIC)
)

##
##
##### Plot of restuls #####
##
##

# re-arranging results table
dfLD <- df00 %>% 
  select(-coef_CS, -coef_CS_low, -coef_CS_high,
         -coef_FN, -coef_FN_low, -coef_FN_high) %>% 
  rename(value = coef_LD, lower = coef_LD_low, upper = coef_LD_high) %>% 
  mutate(coef = "LD")

dfCS <- df00 %>% 
  select(-coef_LD, -coef_LD_low, -coef_LD_high,
         -coef_FN, -coef_FN_low, -coef_FN_high) %>% 
  rename(value = coef_CS, lower = coef_CS_low, upper = coef_CS_high) %>% 
  mutate(coef = "CS")

dfFN <- df00 %>% 
  select(-coef_CS, -coef_CS_low, -coef_CS_high,
         -coef_LD, -coef_LD_low, -coef_LD_high) %>% 
  rename(value = coef_FN, lower = coef_FN_low, upper = coef_FN_high) %>% 
  mutate(coef = "FN")

df01 <- rbind(dfLD, dfCS, dfFN)
df01$coef <- ordered(df01$coef, 
                     levels = c("LD", "CS", "FN"))
head(df01)


## plot with all results (not presented in manuscript)
lnCVR_model_selec_plot <- ggplot(data = df01, 
                                 aes(x = model_ID, y = value, fill = coef)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.75), 
                width = 0,
                color = "black") +
  geom_point(position = position_dodge(width = 0.75),
             shape = 21,
             size = 4) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        axis.text.x = element_text("Arial", size = 10, angle = 90)) +
  scale_fill_manual(values = rev(brewer.pal(n = 3, "Set2"))) +
  scale_color_manual(values = rev(brewer.pal(n = 3, "Set2"))) +
  labs(x = "Model Random Effects Structure", y = "lnCVR estimate") +
  guides(size = FALSE,
         fill = guide_legend(title = "Trait"))

##
##
##### Plot of the best model #####
##
##
base_plot <- orchard_plot_PCL(object = model1, 
                              mod = " ", 
                              est_point_size = 5,
                              alpha = 0.5,
                              cb = FALSE,
                              xlab = "log Coefficient of Variation Ratio (lnCVR)",
                              ylab = "Intercept",
                              transfm = "none") +
  scale_fill_manual(values = rev(brewer.pal(n = 3, "Set2"))) +
  scale_color_manual(values = rev(brewer.pal(n = 3, "Set2"))) +
  scale_y_discrete(labels = c("Laying date", "Clutch size", "# fledlgings"))

cvr_plot <- base_plot + 
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 12),
        axis.text.y = element_text("Arial", size = 12, angle = 45),
        axis.title.y = element_blank(),
        legend.position = "none")

ggsave(filename = "./plots/Figure S6b.jpeg", 
       plot = cvr_plot, 
       device = "jpeg", 
       height = 75, 
       width = 150, 
       units = "mm")

##
## plot no raw data
base_plot <- orchard_plot_PCL_noApples(object = model1, 
                                       mod = " ", 
                                       est_point_size = 5,
                                       alpha = 0.5,
                                       cb = FALSE,
                                       xlab = expression(atop("ln Coefficient of Variation Ratio (lnCVR)",
                                                              ln (CV["Urban"]/CV["Non-urban"]))),
                                       ylab = "Intercept",
                                       transfm = "none") +
  scale_fill_manual(values = rev(brewer.pal(n = 3, "Set2"))) +
  scale_color_manual(values = rev(brewer.pal(n = 3, "Set2"))) +
  scale_y_discrete(labels = c("Laying date", "Clutch size", "# fledlgings"))

cvr_plot_trivar <- base_plot + 
  theme(axis.title = element_text("Arial", size = 10),
        axis.text = element_text("Arial", size = 10),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_x_continuous(limits = c(-0.30,0.30),
                     labels = seq(-0.2,0.2, 0.1),
                     breaks = seq(-0.2,0.2, 0.1))

ggsave(filename = "./plots/Figure 2b.jpeg", 
       plot = cvr_plot_trivar, 
       device = "jpeg", 
       height = 90, 
       width = 90, 
       units = "mm")

##
##
##### Plotting correlations between traits #####
##
##

##
## calculating means and SD per study for plot
##
## calculating means and SD per study for plot
df <- data %>%
  dplyr::select(study_ID, trait, lnCVR, lnRR, species_ID) %>% 
  group_by(study_ID, trait) %>% 
  summarise(mean_rr = mean(lnRR),
            sd_rr = sd(lnRR),
            meancvr = mean(lnCVR),
            sd_cvr = sd(lnCVR)) %>% 
  pivot_wider(id_cols = c(study_ID), 
              names_from = "trait", 
              values_from = c("mean_rr", "sd_rr", "meancvr", "sd_cvr"))
names(df)
names(df) <- c("study_ID", 
               #"species_ID", 
               "mean_rr_Laying_date", "mean_rr_Clutch_size", "mean_rr_Fledglings",
               "sd_rr_Laying_date", "sd_rr_Clutch_size", "sd_rr_Fledglings",
               "meancvr_Laying_date", "meancvr_Clutch_size","meancvr_Fledglings",
               "sd_cvr_Laying_date", "sd_cvr_Clutch_size", "sd_cvr_Fledglings")


# correlation in lnCVR between laying date and clutch size 
lnCVR_LD_CS_cor <- ggplot(data = df, 
                          aes(x = meancvr_Laying_date, 
                              y = meancvr_Clutch_size)) +
  theme_bw() +
  theme(axis.title = element_text(size = 15, family = "Arial"),
        axis.text = element_text(size = 10, family = "Arial"),
        panel.grid = element_blank()) +
  labs(x = "lnCVR Laying date", y = "lnCVR Clutch size") +
  geom_vline(xintercept = 0,
             linetype = 2,
             size = 0.5) +  
  geom_hline(yintercept = 0,
             linetype = 2,
             size = 0.5) +
  geom_errorbar(aes(ymin = meancvr_Clutch_size - sd_cvr_Clutch_size,
                    ymax = meancvr_Clutch_size + sd_cvr_Clutch_size),
                color = "black") + 
  geom_errorbarh(aes(xmin = meancvr_Laying_date - sd_cvr_Laying_date,
                     xmax = meancvr_Laying_date + sd_cvr_Laying_date),
                 color = "black") +
  geom_point(size = 5,
             shape = 21,
             fill = "#74c476",
             color = "black",
             alpha = 1) 

# save plot
ggsave(filename = "./plots/Figure 3d.jpeg", 
       plot = lnCVR_LD_CS_cor,
       height = 100, 
       width = 100, 
       device = "jpeg", 
       units = "mm")

# correlation in lnCVR between laying date and number fledglings
lnCVR_LD_NF_cor <- ggplot(data = df, 
                          aes(x = meancvr_Laying_date, 
                              y = meancvr_Fledglings)) +
  theme_bw() +
  theme(axis.title = element_text(size = 15, family = "Arial"),
        axis.text = element_text(size = 10, family = "Arial"),
        panel.grid = element_blank()) +
  labs(x = "lnCVR Laying date", y = "lnCVR # fledglings") +
  geom_vline(xintercept = 0,
             linetype = 2,
             size = 0.5) +  
  geom_hline(yintercept = 0,
             linetype = 2,
             size = 0.5) +
  geom_errorbar(aes(ymin = meancvr_Fledglings - sd_cvr_Fledglings,
                    ymax = meancvr_Fledglings + sd_cvr_Fledglings),
                color = "black") + 
  geom_errorbarh(aes(xmin = meancvr_Laying_date - sd_cvr_Laying_date,
                     xmax = meancvr_Laying_date + sd_cvr_Laying_date),
                 color = "black") +
  geom_point(size = 5,
             shape = 21,
             fill = "#74c476",
             color = "black",
             alpha = 1) 

# save plot
ggsave(filename = "./plots/Figure 3e.jpeg", 
       plot = lnCVR_LD_NF_cor,
       height = 100, 
       width = 100, 
       device = "jpeg", 
       units = "mm")


# correlation in lnCVR between clutch size and number fledglings
lnCVR_CS_NF_cor <- ggplot(data = df, 
                          aes(x = meancvr_Clutch_size, 
                              y = meancvr_Fledglings)) +
  theme_bw() +
  theme(axis.title = element_text(size = 15, family = "Arial"),
        axis.text = element_text(size = 10, family = "Arial"),
        panel.grid = element_blank()) +
  labs(x = "lnCVR Clutch size", y = "lnCVR # fledglings") +
  geom_vline(xintercept = 0,
             linetype = 2,
             size = 0.5) +  
  geom_hline(yintercept = 0,
             linetype = 2,
             size = 0.5) +
  geom_errorbar(aes(ymin = meancvr_Fledglings - sd_cvr_Fledglings,
                    ymax = meancvr_Fledglings + sd_cvr_Fledglings),
                color = "black") + 
  geom_errorbarh(aes(xmin = meancvr_Clutch_size - sd_cvr_Clutch_size,
                     xmax = meancvr_Clutch_size + sd_cvr_Clutch_size),
                 color = "black") +
  geom_point(size = 5,
             shape = 21,
             fill = "#74c476",
             color = "black",
             alpha = 1) 

# save plot
ggsave(filename = "./plots/Figure 3f.jpeg", 
       plot = lnCVR_CS_NF_cor,
       height = 100, 
       width = 100, 
       device = "jpeg", 
       units = "mm")

