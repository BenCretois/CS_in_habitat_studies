#########################################
# Script run the simulation and analyze #
#########################################
source('Script/simul_ava/simul_functions.R')

# Set the dimensions for the simulation
dim = c(100,100)

# Simulate environmental variables
env <- simul_env(dim)

# Simulate data
data <- simul_data(dim,
                   env,
                   beta0_sp = -7,
                   beta0_cs = 2,
                   beta_forests = 2.5,
                   beta_other_gradient = 3,
                   beta_altitude = -2,
                   beta_roads_sp = 4.5,
                   beta_roads_cs = -6, 
                   beta_urb = -3,
                   beta_nicev = 1,
                   plotdat = TRUE)

################################
# Availability with bias model #
################################

# When no biases the model performs very bad -> because there is no / very few availability points
formula_bias0 <- bias ~ d_roads + d_urb + nice_viewpoints

mod0 <- fit_models(data, formula_bias0)
mod0[[1]] + theme_minimal()

###############################################################
# Try out how robust is the bias model when we drop variables #
###############################################################

formula_bias1 <- bias ~ d_roads + d_urb
formula_bias2 <- bias ~ d_roads
formula_bias3 <- bias ~ d_urb + nice_viewpoints
formula_bias4 <- bias ~ d_urb
formula_bias5 <- bias ~ nice_viewpoints

# We can bootstrap the RMSE to be sure of the improvement
boot = 30
a <- boot_rmse(data, formula_bias0)
b <- boot_rmse(data, formula_bias1)
c <- boot_rmse(data, formula_bias2)
d <- boot_rmse(data, formula_bias3)
e <- boot_rmse(data, formula_bias4)
f <- boot_rmse(data, formula_bias5)

# Make a df of the results and plot!
a_df <- a[[1]]; b_df <- b[[1]]; c_df <- c[[1]]; d_df <- d[[1]]; 
df_diff <- bind_rows(a_df, b_df, c_df, d_df) %>% 
  mutate(Bias_model = c('Full', 'd_roads + d_urb', 'd_roads only', 'd_urb + nice_viewpoints'))

ggplot(df_diff, aes(y=Bias_model, col = Bias_model)) +
  geom_point(aes(x = mean)) +
  geom_linerange(aes(xmin = mean - sd, xmax = mean + sd, y = Bias_model)) +
  theme_bw() +
  geom_vline(xintercept = 0)

# The corrected model always improve the RMSE & the coefficient values
# However, if we forget to include the bias variables which also have ecological
# meaning then this leads to serious bias in this variable.

# IF an environmental variables is correlated to d_roads then it might be 
# a HUGE problem in the non-corrected model

# -> Can try to complicate a bit the simulation

# MOST IMPORTANT VARIABLES IS VARIABLES SHARED BY BOTH CS AND SPECIES ->
# A VARIABLE HAVING BTH AN ECOLOGICAL MEANING AND A BIAS MEANING
