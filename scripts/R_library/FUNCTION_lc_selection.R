##
##
## Function to select land cover data for different buffer sizes
lc_selection <- function(df, buffer_size) {
  
  # creating new variables
  df[["urban_urban_value"]] <- NA
  df[["urban_het"]] <- NA
  df[["rural_urban_value"]] <- NA
  df[["rural_het"]] <- NA
  
  for(i in 1:nrow(df)){
      df[["urban_urban_value"]][i] <- df[[paste0("urban_urbanization_", buffer_size)]][i]
      df[["urban_het"]][i] <- df[[paste0("urban_heterogeneity_", buffer_size)]][i]
      df[["nonurban_urban_value"]][i] <- df[[paste0("nonurban_urbanization_", buffer_size)]][i]
      df[["nonurban_het"]][i] <- df[[paste0("nonurban_heterogeneity_", buffer_size)]][i]
  }
  return(df)
}