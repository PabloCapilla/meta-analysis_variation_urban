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
#' This script compares effect sizes originally extracted by PC-L with those re-extracted by MJT as
#' as part of effect size extraction validation
#' 
##
##
##

##
##### libraries #####
##
pacman::p_load(openxlsx, dplyr, tidyr, metafor) 

##
##### data #####
##

##
## dataset with original and re-extracted effect sizes
data <- read.xlsx("./literature_search/data_extraction_validation/META-ANALYSIS - URBAN_RURAL VARIATION - VALIDATION.xlsx",
          colNames=T,
          sheet = 2)
head(data)
nrow(data)
length(unique(data$study_ID))

##
##
##### Comparison of effect sizes #####
##
##

# re-format dataframe for easy comparisong
data_reformat <- data %>% 
  pivot_longer(cols = 6:ncol(data), names_to = "variable", values_to = "value") %>% 
  separate(col = variable, sep = "_", into = c("type", "habitat", "effect_size")) %>% 
  pivot_wider(names_from = type, values_from = value)

## comparison of means
# plot
ggplot(data = data_reformat %>% 
         filter(effect_size == "mean"), 
       aes(x= original, y = reextracted)) +
  geom_point() +
  theme_bw() +
  stat_smooth(method = "lm")

# correlation
cor.test(data_reformat$original[data_reformat$effect_size == "mean"],
         data_reformat$reextracted[data_reformat$effect_size == "mean"])


## comparison of sd
# plot
ggplot(data = data_reformat %>% 
         filter(effect_size == "sd"), 
       aes(x= original, y = reextracted)) +
  geom_point() +
  theme_bw() +
  stat_smooth(method = "lm")

# correlation
cor.test(data_reformat$original[data_reformat$effect_size == "sd"],
         data_reformat$reextracted[data_reformat$effect_size == "sd"])
  

## comparison of sd without CH0009
# plot
ggplot(data = data_reformat %>% 
         filter(study_ID != "CH0009") %>% 
         filter(effect_size == "sd"), 
       aes(x= original, y = reextracted)) +
  geom_point() +
  theme_bw() +
  stat_smooth(method = "lm")

# correlation
cor.test(data_reformat$original[data_reformat$effect_size == "sd" & data_reformat$study_ID != "CH0009"],
         data_reformat$reextracted[data_reformat$effect_size == "sd" & data_reformat$study_ID != "CH0009"])


## comparison of sd
# plot
ggplot(data = data_reformat %>% 
         filter(effect_size == "n"), 
       aes(x= original, y = reextracted)) +
  geom_point() +
  theme_bw() +
  stat_smooth(method = "lm")

# correlation
cor.test(data_reformat$original[data_reformat$effect_size == "n"],
         data_reformat$reextracted[data_reformat$effect_size == "n"]) 





