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
#' This script prepares the original dataset for meta-analysis.
#' 
##
##
##

##
##### libraries #####
##

pacman::p_load(openxlsx, dplyr, tidyr, metafor, ggplot2, extrafont, RColorBrewer) 
loadfonts()

##
##### data #####
##

##
## study specific data
##
st_df  <- read.xlsx("data/META-ANALYSIS - URBAN_RURAL VARIATION - DATA BASE.xlsx",
                        colNames=T,
                        sheet = 1)
head(st_df)

## initial data exploration

# exploring Country
table(st_df$Country) 
table(is.na(st_df$Country)) 

# exploring study_class
table(st_df$study_class)
table(is.na(st_df$study_class)) 


##
## species specific data 
##
sp_df  <- read.xlsx("data/META-ANALYSIS - URBAN_RURAL VARIATION - DATA BASE.xlsx",
                    colNames=T,
                    sheet = 3)
head(sp_df)

#exploring variables
table(sp_df$scientific_name)
table(sp_df$common_name)
table(sp_df$genus)
table(sp_df$family)
table(sp_df$order)


##
## Comparisons for urban and non-urban populations (effect size data)
##
es_df <- read.xlsx("data/META-ANALYSIS - URBAN_RURAL VARIATION - DATA BASE.xlsx",
                   colNames=T,
                   sheet = 2, )
head(es_df)

#exploring variables
table(es_df$taxa)
table(es_df$species_ID)
table(es_df$scientific_name)
table(es_df$common_name)
table(es_df$trait_class)
table(es_df$trait)
table(es_df$units) # UNKNOWN VALUES OK - CHECKED

##
##### Combining st_df (study-specific data), sp_df (species-specific data) and effect size table #####
##
data01 <- left_join(x = es_df, 
                    y = st_df %>% select(study_ID, 
                                         first_year, 
                                         publication_year, 
                                         urban_lat,
                                         urban_lon,
                                         exact_urban_coor,
                                         rural_lat,
                                         rural_lon,
                                         exact_rural_coor,
                                         study_class, 
                                         Country,
                                         Pop_ID),
                    by = "study_ID")
data <- left_join(x = data01, 
                  y = sp_df %>% select(species_ID, 
                                       family),
                  by = "species_ID")

head(data)

# remove these data columns to have a cleaner data set for analysis
data$details_ES <- NULL


##
##### Exploring effect sizes #####
##
summary(data$urban_mean)
summary(data$urban_sd)
summary(data$urban_se)
summary(data$urban_n)
summary(data$rural_mean)
summary(data$rural_sd)
summary(data$rural_se)
summary(data$rural_n)

table(is.na(data$trait))

# exploring the NA's for the variables of interest
data[is.na(data$urban_mean),]
# OK

data[is.na(data$urban_sd),]
# it seems that 6 papers did not report urban_SD but 2 of those papers did report rural_SD but not urban_SD
## Yes, it is correct, double-checked

data[(is.na(data$urban_se)),]
# OK

data[is.na(data$urban_n),]
# OK

data[is.na(data$rural_mean),]
# OK

data[is.na(data$rural_sd),]
# it seems that 6 papers did not report rural_SD but 2 of those paper did report urban_SD but not rural_SD
## Yes, it is correct (n = 1 nest in both cases) - double-checked

data[(is.na(data$rural_se)),]
# OK

data[is.na(data$rural_n),]
# OK

# double-checking se to sd transformation urban data
round(data[!(is.na(data$urban_se)),"urban_sd"],2) == round((data[!(is.na(data$urban_se)),"urban_se"]*sqrt(data[!(is.na(data$urban_se)),"urban_n"])),2)
# double-checking se to sd transformation non-urban data
round(data[!(is.na(data$urban_se)),"rural_sd"],2) == round((data[!(is.na(data$urban_se)),"rural_se"]*sqrt(data[!(is.na(data$urban_se)),"rural_n"])),2)


# exploring sd/se = 0. 
data[data$urban_sd==0 & !(is.na(data$urban_sd)),] # CHECKED - authors report no variation
data[data$urban_se==0 & !(is.na(data$urban_se)),] # CHECKED ABOVE
data[data$rural_sd==0 & !(is.na(data$rural_sd)),] # CHECKED - authors report no variation
data[data$rural_se==0 & !(is.na(data$rural_se)),] # NONE

##
##
##### Remove observations with n = NA and n = 1 #####
##
##
data <- data %>% 
  filter(!is.na(urban_n)) %>% 
  filter(!is.na(rural_n))

nrow(data[which(data$urban_n == 1 | data$rural_n == 1),]) # 11 obs with n = 1
data <- data[data$urban_n != 1 & data$rural_n != 1,]



##
##
##### Imputing SD values for SD = 0 and for SD = NA #####
##
##

# variable to record imputation method
data$sd_imputation_urban <- "NONE"
data$sd_imputation_rural <- "NONE"


## are there sd/mean values = 0? - # I will remove these further down
data[data$urban_mean==0 & !(is.na(data$urban_mean)),] # no
data[data$rural_mean==0 & !(is.na(data$rural_mean)),] # 1
data[data$rural_sd==0 & !(is.na(data$rural_sd)),] # 2
data[data$urban_sd==0 & !(is.na(data$urban_sd)),] #2

## removing cases of SD = 0 
#data <- data[data$urban_sd != 0 | data$rural_sd != 0,] 

##
## imputation for SD = NA
##

##
## Correlations between mean and SD
##
df_cor <- data.frame(trait = c(data$trait,data$trait),
                     mean_value = c(data$rural_mean,data$urban_mean),
                     sd_value = c(data$rural_sd, data$urban_sd),
                     habitat = c(rep("non_urban", length = length(data$rural_mean)),
                                 rep("urban", length = length(data$urban_mean)))) 

## number of NA to impute per trait
df_cor %>% 
  group_by(trait, habitat) %>% 
  filter(is.na(sd_value)) %>% 
  summarise(n_na = n())

## plots correlations mean-sd and imputation

# Laying date
m_LD_df <- df_cor %>% 
  filter(trait == "laying_date") %>% 
  filter(!is.na(sd_value))

# LD SD data to be predicted
m_LD_pred <- df_cor %>% # data to be predicted
  filter(trait == "laying_date") %>% 
  filter(is.na(sd_value))

#model LD
m_LD <- lm(log(sd_value) ~ log(mean_value),  
           data = m_LD_df %>% 
             filter(sd_value > 0.01))
summary(m_LD)

m_LD_df$predicted <- "NO"
m_LD_pred$predicted <- "YES"
m_LD_pred$sd_value <- exp(predict(object = m_LD, newdata = m_LD_pred %>% select(mean_value)))


plot_LD <- ggplot(data = df_cor %>%
         filter(trait == "laying_date") %>%
         filter(!is.na(sd_value)) %>% 
         group_by(trait) %>% 
         mutate(trait = "Laying date",
                mean_value = log(mean_value),
                sd_value = log(sd_value)), 
       aes(x = mean_value, y = sd_value, fill = habitat)) +
  #facet_grid(~habitat)+
  geom_point(color = "black", size = 2, shape = 21,
             fill = brewer.pal(n = 3, name = "Set2")[3]) +
  facet_wrap(~trait) +
  stat_smooth(method = "lm", formula = y~x, color = "black", fill = "grey") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text("Arial", size = 10),
        axis.text = element_text("Arial", size = 10),
        strip.background = element_blank(),
        strip.text = element_text("Arial", size = 10)) +
  labs(x = expression(ln(bar(x))), y = "ln(SD)") +
  scale_x_continuous(breaks = seq(3,6,0.5),
                     labels = seq(3,6,0.5)) +
  scale_y_continuous(breaks = seq(0,4,1),
                     labels = seq(0,4,1),
                     limits = c(-0.5,4)) 

plot_LD_pred <- plot_LD +
  geom_point(data = m_LD_pred %>% 
               mutate(trait = "Laying date",
                      mean_value = log(mean_value),
                      sd_value = log(sd_value)),
             fill = brewer.pal(n = 3, name = "Set2")[3],
             color = "black", shape = 23, size = 2.5)


ggsave(filename = "./plots/Figure S2a.jpeg",
       plot = plot_LD_pred, 
       device = "jpeg",
       units = "mm", 
       width = 90, 
       height = 90)


## predictions to data base
data$sd_imputation_urban[is.na(data$urban_sd) & data$trait == "laying_date"] <- "lm"
data$urban_sd[is.na(data$urban_sd) & data$trait == "laying_date"] <- exp(predict(object = m_LD, 
                                                                                 newdata = data.frame(mean_value = data$urban_mean[is.na(data$urban_sd) & data$trait == "laying_date"])))

data$sd_imputation_rural[is.na(data$rural_sd) & data$trait == "laying_date"] <- "lm"
data$rural_sd[is.na(data$rural_sd) & data$trait == "laying_date"] <- exp(predict(object = m_LD, 
                                                                                 newdata = data.frame(mean_value = data$rural_mean[is.na(data$rural_sd) & data$trait == "laying_date"])))

table(data$sd_imputation_rural[data$trait == "laying_date"], 
      data$study_ID[data$trait == "laying_date"])
table(data$sd_imputation_urban[data$trait == "laying_date"],
      data$study_ID[data$trait == "laying_date"])


##
## Clutch size (CS)

# data CS
table(is.na(df_cor$sd_value))


m_CS_df <- df_cor %>% 
  filter(trait == "clutch_size") %>% 
  filter(!is.na(sd_value))

m_CS_pred <- df_cor %>% # data to be predicted
  filter(trait == "clutch_size") %>% 
  filter(is.na(sd_value))

# model CS
m_CS <- lm(log(sd_value) ~ log(mean_value),  
           data = m_CS_df %>% 
             filter(sd_value > 0.01))
summary(m_CS)

m_CS_df$predicted <- "NO"
m_CS_pred$predicted <- "YES"
m_CS_pred$sd_value <- exp(predict(object = m_CS, newdata = m_CS_pred %>% select(mean_value)))

plot_CS <- ggplot(data = df_cor %>%
                    filter(trait == "clutch_size") %>%
                    filter(!is.na(sd_value)) %>% 
                    group_by(trait) %>% 
                    mutate(trait = "Clutch size",
                           mean_value = log(mean_value),
                           sd_value = log(sd_value)), 
                  aes(x = mean_value, y = sd_value)) +
  geom_point(fill = brewer.pal(n = 3, name = "Set2")[2], 
              color = "black", size = 2, shape = 21) +
  facet_wrap(~trait) +
  stat_smooth(data = df_cor %>%
                filter(trait == "clutch_size") %>%
                filter(!is.na(sd_value)) %>% 
                group_by(trait) %>% 
                mutate(trait = "Clutch size",
                       mean_value = log(mean_value),
                       sd_value = log(sd_value)), 
              method = "lm", formula = y~x, color = "black", fill = "grey") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text("Arial", size = 10),
        axis.text = element_text("Arial", size = 10),
        strip.background = element_blank(),
        strip.text = element_text("Arial", size = 10)) +
  labs(x = expression(ln(bar(x))), y = "ln(SD)") +
  scale_x_continuous(breaks = seq(0.5,3,0.5),
                     labels = seq(0.5,3,0.5)) +
  scale_y_continuous(breaks = seq(-1,2,1),
                     labels = seq(-1,2,1),
                     limits = c(-1,2)) 
plot_CS_pred <- plot_CS +
  geom_point(data = m_CS_pred %>% 
               mutate(trait = "Clutch size",
                      trait = "Clutch size",
                      mean_value = log(mean_value),
                      sd_value = log(sd_value)),
             fill = brewer.pal(n = 3, name = "Set2")[2],
             color = "black", shape = 23, size = 2.5)

ggsave(filename = "./plots/Figure S2b.jpeg",
       plot = plot_CS_pred, 
       device = "jpeg",
       units = "mm", 
       width = 90, 
       height = 90)


## predictions to data base
data$sd_imputation_urban[is.na(data$urban_sd) & data$trait == "clutch_size"] <- "lm"
data$urban_sd[is.na(data$urban_sd) & data$trait == "clutch_size"] <- exp(predict(object = m_CS, 
                                                                                 newdata = data.frame(mean_value = data$urban_mean[is.na(data$urban_sd) & data$trait == "clutch_size"])))

data$sd_imputation_rural[is.na(data$rural_sd) & data$trait == "clutch_size"] <- "lm"
data$rural_sd[is.na(data$rural_sd) & data$trait == "clutch_size"] <- exp(predict(object = m_CS, 
                                                                                 newdata = data.frame(mean_value = data$rural_mean[is.na(data$rural_sd) & data$trait == "clutch_size"])))
table(data$sd_imputation_rural)
table(data$sd_imputation_urban)

table(data$sd_imputation_rural[data$trait == "clutch_size"], 
      data$study_ID[data$trait == "clutch_size"])
table(data$sd_imputation_urban[data$trait == "clutch_size"],
      data$study_ID[data$trait == "clutch_size"])

# Number of fledglings (FN)
# data FN
table(is.na(df_cor$sd_value))

m_FN_df <- df_cor %>% 
  filter(trait == "fledglings_per_attempt") %>% 
  filter(!is.na(sd_value))

m_FN_pred <- df_cor %>% # data to be predicted
  filter(trait == "fledglings_per_attempt") %>% 
  filter(is.na(sd_value))

# model FN
m_FN <- lm(log(sd_value) ~ log(mean_value),  
           data = m_FN_df %>% 
             filter(sd_value > 0.01))
summary(m_FN)

m_FN_df$predicted <- "NO"
m_FN_pred$predicted <- "YES"
m_FN_pred$sd_value <- exp(predict(object = m_FN, newdata = m_FN_pred %>% select(mean_value)))

plot_FN <- ggplot(data = m_FN_df %>%
                    filter(sd_value != 0.01) %>% 
                    group_by(trait) %>% 
                    mutate(trait = "Number of fledglings",
                           mean_value = log(mean_value),
                           sd_value = log(sd_value)), 
                  aes(x = mean_value, y = sd_value)) +
  geom_point(fill = brewer.pal(n = 3, name = "Set2")[1], 
             color = "black", size = 2, shape = 21) +
  facet_wrap(~trait) +
  stat_smooth(data = m_FN_df %>%
                filter(sd_value != 0.01) %>% 
                group_by(trait) %>% 
                mutate(trait = "Number of fledglings",
                       mean_value = log(mean_value),
                       sd_value = log(sd_value)), 
              method = "lm", formula = y~x, color = "black", fill = "grey") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text("Arial", size = 10),
        axis.text = element_text("Arial", size = 10),
        strip.background = element_blank(),
        strip.text = element_text("Arial", size = 10)) +
  labs(x = expression(ln(bar(x))), y = "ln(SD)") +
  scale_x_continuous(breaks = seq(-2,3,1),
                     labels = seq(-2,3,1)) +
  scale_y_continuous(breaks = seq(-1,2,1),
                     labels = seq(-1,2,1),
                     limits = c(-1.5,2)) 

plot_FN_pred <- plot_FN +
  geom_point(data = m_FN_pred %>% 
               mutate(trait = "Number of fledglings",
                      mean_value = log(mean_value),
                      sd_value = log(sd_value)),
             fill = brewer.pal(n = 3, name = "Set2")[1], 
             color = "black", shape = 23, size = 2.5)

ggsave(filename = "./plots/Figure S2c.jpeg",
       plot = plot_FN_pred, 
       device = "jpeg",
       units = "mm", 
       width = 90, 
       height = 90)




## predictions to data base
data$sd_imputation_urban[is.na(data$urban_sd) & data$trait == "fledglings_per_attempt"] <- "lm"
data$urban_sd[is.na(data$urban_sd) & data$trait == "fledglings_per_attempt"] <- exp(predict(object = m_FN, 
                                                                                 newdata = data.frame(mean_value = data$urban_mean[is.na(data$urban_sd) & data$trait == "fledglings_per_attempt"])))

data$sd_imputation_rural[is.na(data$rural_sd) & data$trait == "fledglings_per_attempt"] <- "lm"
data$rural_sd[is.na(data$rural_sd) & data$trait == "fledglings_per_attempt"] <- exp(predict(object = m_FN, 
                                                                                 newdata = data.frame(mean_value = data$rural_mean[is.na(data$rural_sd) & data$trait == "fledglings_per_attempt"])))

table(data$sd_imputation_rural[data$trait == "fledglings_per_attempt"])
table(data$sd_imputation_urban[data$trait == "fledglings_per_attempt"])

table(data$sd_imputation_rural[data$trait == "fledglings_per_attempt"], 
      data$study_ID[data$trait == "fledglings_per_attempt"])
table(data$sd_imputation_urban[data$trait == "fledglings_per_attempt"],
      data$study_ID[data$trait == "fledglings_per_attempt"])

table(data$sd_imputation_rural, data$study_ID)
table(data$sd_imputation_urban, data$study_ID)

##
##### Descriptive summary of data (whole data set) #####
##

# removing obs with sd = 0
data <- data %>% 
  filter(!(urban_sd == 0)) %>% 
  filter(!(rural_sd == 0)) 
         

# Missing data prenst in means, SD or N? (all should have been removed already)
table(is.na(data$urban_mean))
table(is.na(data$urban_sd))
table(is.na(data$urban_n)) 
table(is.na(data$rural_mean))
table(is.na(data$rural_sd))
table(is.na(data$rural_n))

# no observations with sd = 0
data[data$urban_sd ==0 | data$rural_sd ==0,]
 
##
##
##### Summary of final data set #####
##
##

nrow(data) # number of ES
length(unique(data$species_ID)) # number of species
length(unique(data$study_ID))   # number of studies included


# number of studies with ES per trait
data %>% 
  group_by(trait) %>% 
  summarise(n_cases_per_study = n()) 

# number of studies with ES per trait
data %>% 
  group_by(trait, study_ID) %>% 
  summarise(n_cases_per_study = n()) %>% 
  group_by(trait) %>% 
  summarise(n_studies = n())

# number of ES per study
class(data$trait)

# single and multiple year studies
data %>% 
  group_by(multiple_year_ES, study_ID) %>% 
  summarise(n_cases_per_study = n()) %>% 
  group_by(multiple_year_ES) %>% 
  summarise(n_studies = n())

nrow(data %>% filter(multiple_year_ES == "YES"))
nrow(data %>% filter(multiple_year_ES == "NO"))

# observations for within and among-year analysis
table(data$multiple_year_ES)

data$trait <- as.factor(data$trait)
table(data$trait)
levels(data$trait) <- c("Clutch size", "# Fledglings", "Laying date")
table(data$trait)

data$trait <- ordered(data$trait, 
                      levels = c("Laying date", "Clutch size", "# Fledglings"))
table(data$trait)

data %>% 
  group_by(trait, study_ID) %>% 
  summarise(n_cases_per_study = n()) %>% 
  summarise(sum(n_cases_per_study))

obs_study <- data %>% 
  group_by(trait, study_ID) %>% 
  summarise(n_cases_per_study = n()) %>% 
  ggplot(aes(x = n_cases_per_study, fill = trait)) +
  scale_fill_manual(values = rev(brewer.pal(n = 3, name = "Set2"))) +
  geom_histogram(color = "black") +
  facet_grid(~trait) +
  theme_bw() +
  labs(x = "Number of urban - non-urban comparisons", y = "Study count") +
  theme(strip.text = element_text("Arial", size = 10),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title = element_text("Arial", size = 10),
        axis.text = element_text("Arial", size = 8),
        panel.grid = element_blank()) +
  scale_y_continuous(breaks = seq(0,25,5),
                     labels = seq(0,25,5),
                     limits = c(0,25)) 

ggsave(filename = "./plots/Figure S3.jpeg",
       plot = obs_study, 
       device = "jpeg",
       units = "mm", 
       width = 120, 
       height = 90)
##
##### Calculating meta-analysis y variables #####
##

## lnCVR
data <- as.data.frame(escalc(measure="CVR", 
                             n1i=urban_n, n2i=rural_n, 
                             m1i=urban_mean, m2i=rural_mean, 
                             sd1i=urban_sd, sd2i=rural_sd,
                             data=data,
                             var.names=c("lnCVR","lnCVR.sv"), 
                             add.measure=F,
                             append=TRUE))

## lnVR
data <- as.data.frame(escalc(measure="VR", 
                             n1i=urban_n, n2i=rural_n, 
                             m1i=urban_mean, m2i=rural_mean, 
                             sd1i=urban_sd, sd2i=rural_sd,
                             data=data,
                             var.names=c("lnVR","lnVR.sv"), 
                             add.measure=F,
                             append=TRUE))

## lnRR
data <- as.data.frame(escalc(measure="ROM", 
                             n1i=urban_n, n2i=rural_n, 
                             m1i=urban_mean, m2i=rural_mean, 
                             sd1i=urban_sd, sd2i=rural_sd,
                             data=data,
                             var.names=c("lnRR","lnRR.sv"), 
                             add.measure=F,
                             append=TRUE))
head(data)



##
##
##### Saving full table for analysis #####
saveRDS(object = data, file = "./data/processed_RDS_data_files/metaanalysis_202220318.RDS")

