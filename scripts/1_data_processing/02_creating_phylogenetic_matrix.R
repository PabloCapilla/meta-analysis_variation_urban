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
#' This script creates a matrix of phylogenetic correlation to include in meta-analytical models
#' 
##
##
##


##
##### libraries #####
##
pacman::p_load(dplyr, tidyr, metafor, rotl, 
               ape, curl, treeio, ggtree,
               phytools, RColorBrewer, ggimage, extrafont)
loadfonts()

##
##### data #####
##
data <- readRDS("./data/processed_RDS_data_files/metaanalysis_full_data.RDS")
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
##### Tree plot and proportions of data per species #####
##
##
phylo_bar <- data %>% 
  group_by(scientific_name, trait, family) %>% 
  summarise(n_obs = n()) %>% 
  spread(key = "trait", value = "n_obs") %>% 
  mutate(total_obs = sum(`Laying date`, `Clutch size`, `# Fledglings`, na.rm = T)) %>% 
  mutate(`Laying date` = `Laying date` / total_obs,
         `Clutch size` = `Clutch size` / total_obs,
         `# Fledglings` = `# Fledglings` / total_obs,
         scientific_name =sub(x = scientific_name, pattern = " ", replacement = "_"))
phylo_bar$scientific_name <- recode(phylo_bar$scientific_name,
                                    "Accipiter_cooperi" = "Accipiter_cooperii",
                                    "Athene_cunicuaria" = "Athene_cunicularia")

# sort and re-arrange dataframe
phylo_bar[is.na(phylo_bar)] <- 0
phylo_bar_plot <- phylo_bar %>% 
  select(scientific_name, `Laying date`, `Clutch size`, `# Fledglings`)

phylo_bar <- as.data.frame(phylo_bar)
phylo_bar_plot <- as.data.frame(phylo_bar_plot)
row.names(phylo_bar_plot) <- phylo_bar_plot$scientific_name
phylo_bar_plot$scientific_name <- NULL

phylo_bar_plot2 <- phylo_bar_plot %>% 
  mutate(label = row.names(phylo_bar_plot)) %>% 
  gather(key = "trait", value = "value", 1:3)
phylo_bar_plot2$family <- NULL
phylo_bar_plot2$trait <- ordered(phylo_bar_plot2$trait, 
                                 levels = c("# Fledglings", "Clutch size", "Laying date"))

##
## phylogenetic tree with sample size per species
p2 <- ggtree(tree, size = 1.2) %<+% phylo_bar + 
  xlim_tree(c(0, 35)) +
  geom_tiplab(aes(label = paste0("italic('", sub(x = label, pattern = "_", replacement = " "), "')")), 
              parse = TRUE, 
              offset = 0.5,
              size = 2.8,
              hjust = 0) +
  geom_text(aes(label = total_obs, x = 36), size = 2.5)
 
##
## add panel with proportion of data per species and trait
panel <- facet_plot(p2, 
           panel="Proportion of observations", 
           data=phylo_bar_plot2,
           ggstance::geom_colh, 
           mapping = aes(x = value, fill = trait, color = trait),
           position="stackv") +
  theme(legend.position = "bottom") +
  theme(strip.text = element_text(size = 10, "Arial"),
        strip.background = element_blank(),
        legend.text = element_text(size = 7, "Arial"),
        legend.title = element_blank()) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Set2")) +
  scale_color_manual(values = brewer.pal(n = 3, name = "Set2"))


ggsave(filename = "./plots/Figure 1a.jpeg",
       height = 110, 
       plot = panel, 
       device = "jpeg", 
       width = 150, 
       units = "mm")


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
saveRDS(phylo_cor, file = "./data/processed_RDS_data_files/phylogenetic_correlations.RDS")



