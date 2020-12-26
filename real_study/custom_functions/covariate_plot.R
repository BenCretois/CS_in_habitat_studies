################################################################################
# Script to plot the covariate distribution for CS, telemetry and availability #
################################################################################

# The function returns a boxplot of the distribution of the covariates for each 
# "type" of observations: Telemetry, Citizen Science and availability

covariate_plot <- function(data){
  
  # Take only telem, CS and availability
  data <- data %>% 
    filter(type == 'telem_gps' | type == 'cs' | type == 'random_cs')
  
  # Make the plot
  data_descriptive <- data %>% st_drop_geometry() %>%  
    dplyr::select(type, path_use_log, d_roads, d_urb, n_forest, n_agr, slope, altitude) %>% 
    gather(key = 'predictors', value = 'predictors_value', -type)
  
  # plot boxplots of all variables for all "types"
  cov_plots <- ggplot(data_descriptive, aes(x = type, y = predictors_value, fill = type)) +
    geom_boxplot() +
    geom_smooth() + 
    theme_bw() + 
    facet_wrap(~predictors, scales = "free_y") +
    scale_fill_viridis(option = 'D', discrete = TRUE, alpha = .7)
  
  # Object to return
  return(cov_plots)
}
