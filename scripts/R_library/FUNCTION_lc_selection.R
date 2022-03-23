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
      df[["urban_urban_value"]][i] <- df[[paste0("urbanization.", buffer_size,"_urban")]][i]
      df[["urban_het"]][i] <- df[[paste0("heterogeneity.", buffer_size,"_urban")]][i]
      df[["rural_urban_value"]][i] <- df[[paste0("urbanization.", buffer_size,"_rural")]][i]
      df[["rural_het"]][i] <- df[[paste0("heterogeneity.", buffer_size,"_rural")]][i]
  }
  return(df)
}