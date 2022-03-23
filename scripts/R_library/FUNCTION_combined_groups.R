##
## function to combined means and sds across groups of observations
##
combine_groups <- function(means, sds, ns) {
  
  ## compute combined means
  combined_mean <- sum((means*ns))/sum(ns)
  
  ## compute combined sd
  dev_means <- sum(((means - combined_mean)^2)*ns) # deviation from group means
  dev_sd <- sum(sds^2*(ns)) # sum of group sd
  combined_sd <- sqrt((dev_means + dev_sd)/(sum(ns))) # combined sd
  
  return(c(combined_mean = combined_mean, combined_sd = combined_sd))
}
