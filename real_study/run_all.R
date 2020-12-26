###############################################################################
################ Script for the analysis CS vs Telemetry data #################
###############################################################################

####################
# SOME INFORMATION #
####################

# For all species I selected telemetry observations and Citizen Science observations during 
# ONLY summer and ONLY during the day!

###################
# Set up analysis #
###################

# Load the required function and libraries
source('analysis/real_analysis/scripts/Libraries.R')

# Source all the functions
sourceFolder <- function(folder, recursive = FALSE, ...) 
{ 
  files <- list.files(folder, pattern = "[.][rR]$", 
                      full.names = TRUE, recursive = recursive)
  if (!length(files))
    stop(simpleError(sprintf('No R files in folder "%s"', folder)))
  src <- invisible(lapply(files, source, ...))
  message(sprintf('%s files sourced from folder "%s"', length(src), folder))
}

sourceFolder('analysis/real_analysis/scripts/custom_functions', recursive = TRUE)


tmap_mode('view')
set.seed(7)

####################
# Make the dataset #
####################

# Get a standardized dataset for all 3 species. Thinning reduce the amount of neighbouring telemetry data
# reducing the risk of any temporal / spatial correlation.

moose <- make_moose_dataset(buffer_for_cs = 10000, buffer_for_telem = 5000, thinning = hours(5), ratio_random_cs = 5)
roedeer <- make_roedeer_dataset(buffer_for_cs = 10000, buffer_for_telem = 5000, thinning = hours(5), ratio_random_cs = 5)
wreindeer <- make_reindeer_dataset(buffer_for_cs = 1, buffer_for_telem = 5000, thinning = hours(10), ratio_random_cs = 5)

# Check how many observations
table(moose$type)
table(roedeer$type)
table(wreindeer$type)

######################################
# Locations of CS and Telemetry data #
######################################
moose_map <- moose %>% filter(type == 'telem_gps' | type == 'cs')
tm_shape(moose_map) + tm_dots(col = 'type')

roedeer_map <- roedeer %>% filter(type == 'telem_gps' | type == 'cs')
tm_shape(roedeer_map) + tm_dots(col = 'type')

wreindeer_map <- wreindeer %>% filter(type == 'telem_gps' | type == 'cs')
tm_shape(wreindeer_map) + tm_dots(col = 'type')

#############################################
# Get the covariate values for each dataset #
#############################################

# List the covariate raster -> call them inside the function
list_raster_mr <- list.files("analysis/real_analysis/covariates_moose_roe/")
list_raster_wreindeer <- list.files("analysis/real_analysis/covariates_reindeer/")

moose <- get_covariate('moose', moose, list_raster_mr, 4)
roedeer <- get_covariate('roedeer', roedeer, list_raster_mr, 4)
wreindeer <- get_covariate('wreindeer', wreindeer, list_raster_wreindeer, 4)

#####################################################################
# Distribution of the covariates for telemetry, CS and availability #
#####################################################################

covariate_plot(moose)
covariate_plot(roedeer)
covariate_plot(wreindeer)

################################################################
# Get the corrected availability points & check the bias model #
################################################################

# Observer model -> moose and roe deer
formula_obs_mr <- type_recode ~  d_roads_log + d_urb_log + path_use_log + n_forest_log + pop_log
formula_obs_reindeer <- type_recode ~ d_roads_log + d_urb_log + path_use_log + n_forest_log 

# Get the corrected availability points
ava_moose <- get_corrected_availability('moose', moose, formula_obs_mr, ratio_availability = 5, list_raster_mr)
ava_moose_points <- ava_moose[[1]] %>% get_covariate('moose', ., list_raster_mr, 3) # Take the availability points and add the covariates

ava_roedeer <- get_corrected_availability('roedeer', roedeer, formula_obs_mr, ratio_availability = 5, list_raster_mr)
ava_roedeer_points <- ava_roedeer[[1]] %>% get_covariate('roedeer', ., list_raster_mr, 3) # Take the availability points and add the covariates

ava_wreindeer <- get_corrected_availability('wreindeer', wreindeer, formula_obs_reindeer, ratio_availability = 5, list_raster_wreindeer)
ava_wreindeer_points <- ava_wreindeer[[1]] %>% get_covariate('wreindeer', ., list_raster_wreindeer, 3) # Take the availability points and add the covariates

# Observer model for moose
Efxplot(ava_moose[[2]]) + theme_bw() + coord_flip(ylim = c(-5,5))

# Observer model for for roe deer
Efxplot(ava_roedeer[[2]]) + theme_bw() + coord_flip(ylim = c(-5,5))

# Observer model for for wild reindeer
Efxplot(ava_wreindeer[[2]]) + theme_bw() + coord_flip(ylim = c(-5,5))

# All observer models together
myINLAplot(list(ava_moose[[2]], ava_roedeer[[2]], ava_wreindeer[[2]]),
           ModelNames = c('Moose', 'Roe deer', 'Wild reindeer'), Intercept = FALSE) + theme_bw() + coord_flip(ylim = c(-1, 1))


######################################################
# Sample availability points from another bias model #
######################################################

# For moose with roe deer and wild reindeer bias model
ava_moose_with_rd_bias <- get_availability_cross_species('moose', moose, 'roedeer', ava_roedeer[[2]], ratio_availability = 5, 
                                                         list_raster = list_raster_mr, formula_obs_mr, type = 'corrected_1') %>% 
  get_covariate('moose', ., list_raster_mr, 3)

ava_moose_with_wr_bias <- get_availability_cross_species(species = 'moose', data = moose, species_bias_model = 'wreindeer', ava_wreindeer[[2]], ratio_availability = 5, 
                                                         list_raster = list_raster_mr, formula_obs_reindeer, type = 'corrected_2') %>% 
  get_covariate('moose', ., list_raster_mr, 3)

# For roe deer with moose and wild reindeer bias model
ava_rd_with_moose_bias <- get_availability_cross_species('roedeer', roedeer, 'moose', ava_moose[[2]], ratio_availability = 5, 
                                                         list_raster = list_raster_mr, formula_obs_mr, type = 'corrected_1') %>% 
  get_covariate('roedeer', ., list_raster_mr, 3)

ava_rd_with_wr_bias <- get_availability_cross_species('roedeer', roedeer, 'wreindeer', ava_wreindeer[[2]], ratio_availability = 5, 
                                                         list_raster = list_raster_mr, formula_obs_reindeer, type = 'corrected_2') %>% 
  get_covariate('roedeer', ., list_raster_mr, 3)

# For wild reindeer with moose and roe deer bias model
ava_wr_with_moose_bias <- get_availability_cross_species('wreindeer', wreindeer, 'moose', ava_moose[[2]], ratio_availability = 5, 
                                                         list_raster = list_raster_wreindeer, formula_obs_mr, type = 'corrected_1') %>% 
  get_covariate('wreindeer', ., list_raster_wreindeer, 3)

ava_wr_with_rd_bias <- get_availability_cross_species('wreindeer', wreindeer, 'roedeer', ava_roedeer[[2]], ratio_availability = 5, 
                                                         list_raster = list_raster_wreindeer, formula_obs_mr, type = 'corrected_2') %>% 
  get_covariate('wreindeer', ., list_raster_wreindeer, 3)


#########################################
# Bind the rows to get the full dataset #
#########################################

moose_corrected <- bind_rows(moose, ava_moose_points,ava_moose_with_rd_bias, ava_moose_with_wr_bias) %>% st_transform(., crs = 25833) %>% drop_na(type)
roedeer_corrected <- bind_rows(roedeer, ava_roedeer_points, ava_rd_with_moose_bias, ava_rd_with_wr_bias) %>% st_transform(., crs = 25833) %>% drop_na(type)
wreindeer_corrected <- bind_rows(wreindeer, ava_wreindeer_points, ava_wr_with_moose_bias, ava_wr_with_rd_bias) %>% st_transform(., crs = 25833) %>% drop_na(type)

##################################################
# Visualize the location of the different points #
##################################################
tm_shape(moose_corrected) + tm_dots(col = 'type') #+ tm_shape(ava_moose_points) + tm_dots() 
tm_shape(roedeer) + tm_dots(col = 'type') #+ tm_shape(ava_roedeer_points) + tm_dots() 
tm_shape(wreindeer) + tm_dots(col = 'type') #+ tm_shape(ava_wreindeer_points) + tm_dots() 

##################
# Fit the models #
##################

mod_moose <- fit_rsf_models(moose_corrected)
mod_roedeer <- fit_rsf_models(roedeer_corrected)
mod_wreindeer <- fit_rsf_models(wreindeer_corrected)

# Save the models
saveRDS(mod_moose, 'analysis/real_analysis/saved_models/mod_moose.rds')
saveRDS(mod_roedeer, 'analysis/real_analysis/saved_models/mod_roedeer.rds')
saveRDS(mod_wreindeer, 'analysis/real_analysis/saved_models/mod_wreindeer.rds')

# Read the models
mod_moose <- readRDS('analysis/real_analysis/saved_models/mod_moose.rds')
mod_roedeer <- readRDS('analysis/real_analysis/saved_models/mod_roedeer.rds')
mod_wreindeer <- readRDS('analysis/real_analysis/saved_models/mod_wreindeer.rds')

# Plot the coefficients
mod_moose[[1]] + theme_bw() + coord_flip(ylim = c(-2,2))
mod_roedeer[[2]] + theme_bw() + coord_flip(ylim = c(-2,2))
mod_wreindeer[[2]] + theme_bw() + coord_flip(ylim = c(-3,3))

myINLAplot(list(mod_moose[[3]], mod_moose[[4]], mod_moose[[5]]), Intercept = FALSE) + theme_bw() #+ guides(colour = FALSE)
myINLAplot(list(mod_roedeer[[3]], mod_roedeer[[4]], mod_roedeer[[5]]), Intercept = FALSE) + theme_bw() #+ guides(colour = FALSE)
myINLAplot(list(mod_wreindeer[[3]], mod_wreindeer[[4]], mod_wreindeer[[5]]), Intercept = FALSE, ModelNames = c('Telemetry', 'Naive CS', 'Corrected CS')) + theme_bw()

cowplot::plot_grid(a,b,c, nrow = 1)
?plot_grid
#########################
# Make suitability maps #
#########################

suit_moose <- make_suitability_map(species = 'moose', data = moose, models = mod_moose, cell_size = 1000)
suit_moose[[3]]
suit_moose[[4]]

suit_roedeer <- make_suitability_map(species = 'roedeer', data = roedeer, model = mod_roedeer, cell_size = 1000)
suit_roedeer[[3]]
suit_roedeer[[4]]

suit_wreindeer <- make_suitability_map(species ='wreindeer', data = wreindeer, model = mod_wreindeer, cell_size = 2000)
suit_wreindeer[[3]]
suit_wreindeer[[4]]


############################################## END OF THE ANALYSIS ###########################################
