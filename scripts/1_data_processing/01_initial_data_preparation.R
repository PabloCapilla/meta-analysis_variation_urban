###
###
#' 
#' Script for:
#' A global meta-analysis reveals higher phenological variation in urban birds than in their non-urban neighbours
#' Capilla-Lasheras et al. 
#' Preprint: https://doi.org/10.1101/2021.09.24.461498
#' 
#' Latest update: 2022/06/28
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
##### Initial sample size #####
##
nrow(data) ## n effect sizes
length(unique(data$species_ID)) ## n species
length(unique(data$study_ID)) ## n studies


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

# removing observations with missing sample size
data <- data %>% 
  filter(!is.na(rural_n) | !is.na(urban_n)) # 440 observations retained (3 observations filtered out)

data <- data %>%           # 429 observations retained (11 observations filtered out)
  filter(urban_n != 1) %>% 
  filter(rural_n != 1)

# exploring the NA's for the variables of interest
data[is.na(data$urban_mean),]
# OK

data[is.na(data$urban_sd),]
data[is.na(data$urban_sd),c("study_ID", "urban_sd", "rural_sd")]
# it seems that 10 papers did not report urban_SD but 2 of those papers did report rural_SD but not urban_SD
## Yes, it is correct, double-checked

data[is.na(data$urban_n),]
# OK

data[is.na(data$rural_mean),]
# OK

data[is.na(data$rural_sd),]
data[is.na(data$rural_sd),c("study_ID", "urban_sd", "rural_sd")]
# it seems that 10 papers did not report rural_SD but 2 of those paper did report urban_SD but not rural_SD
## Yes, it is correct (n = 1 nest in both cases) - double-checked

data[is.na(data$rural_n),]
# OK

# double-checking se to sd transformation urban data
round(data[!(is.na(data$urban_se)),"urban_sd"],2) == round((data[!(is.na(data$urban_se)),"urban_se"]*sqrt(data[!(is.na(data$urban_se)),"urban_n"])),2)
# double-checking se to sd transformation non-urban data
round(data[!(is.na(data$urban_se)),"rural_sd"],2) == round((data[!(is.na(data$urban_se)),"rural_se"]*sqrt(data[!(is.na(data$urban_se)),"rural_n"])),2)

# exploring sd = 0. 
data[data$urban_sd==0 & !(is.na(data$urban_sd)),] # CHECKED - authors report no variation
data[data$rural_sd==0 & !(is.na(data$rural_sd)),] # CHECKED - authors report no variation

# how many observations with SD = 0?
nrow(data %>%                # 4 observations with SD = 0 (all four observations double check by going to original source - see above)
       filter(urban_sd == 0 | rural_sd == 0))

# removing observations with SD = 0
data <- data[-which(data$urban_sd == 0 | data$rural_sd == 0),]  # 425 observations retained (4 observations filtered out)

# how many observations with SD = NA
nrow(data %>%                # 26 observations with SD = NA
       filter(is.na(urban_sd)| is.na(rural_sd)))
data <- data[-which(is.na(data$urban_sd) | is.na(data$rural_sd)),]  # 399 observations retained (26 observations filtered out)


##
##### Descriptive summary of final data #####
##

## description of final data set
##
nrow(data) # number of ES
length(unique(data$species_ID)) # number of species
length(unique(data$study_ID))   # number of studies included

summary(data$urban_mean)
summary(data$urban_sd)
summary(data$urban_n)
summary(data$rural_mean)
summary(data$rural_sd)
summary(data$rural_n)


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

ggsave(filename = "./plots/Figure S2.jpeg",
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

## SMDH
data <- as.data.frame(escalc(measure="SMDH", 
                             n1i=urban_n, n2i=rural_n, 
                             m1i=urban_mean, m2i=rural_mean, 
                             sd1i=urban_sd, sd2i=rural_sd,
                             data=data,
                             var.names=c("SMDH","SMDH.sv"), 
                             add.measure=F,
                             append=TRUE))



head(data)

##
##### Save full table #####
saveRDS(object = data, file = "./data/processed_RDS_data_files/metaanalysis_full_data.RDS")
#####













