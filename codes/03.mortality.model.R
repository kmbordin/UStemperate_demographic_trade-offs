########################################################################
####### MORTALITY ESTIMATES ACROSS SPECIES #############################
########################################################################

# Mortality data -----  
# DO NOT RUN  - lazy and large memory requirement; use the RData! -----
# load("rdata/live.RData")
# mort <- live #good data
# data_200 <- species_list_to_mortality(x = mort,min.n.surv = 200, min.n.mort = 5)
# fit.stan <- data_200

# Load the data
# load("rdata/fit.stan.total.rdata")
# 
# Growth is already transformed. 
# Re-generate the raw annualised absolute growth rate values:
# 
# data <- fit.stan
# str(data)
# data$gr_raw <- data$gr ^ (1/0.47)
# 
# Explore whether there are potential positive growth outliers (errors / impossible / extremely unlikely growth values):
# 
# data %>% 
#   ggplot(aes(exp(dbh), gr_raw)) +
#   geom_point() +
#   facet_wrap(~sp) +
#   labs(x = "dbh (cm)",
#        y = "growth (cm yr-1)") +
#   mytheme
# 
# data2 = data %>% filter(!is.na(gr_raw)) 
# 
# # Quality checks  = removing outliers
# data2 <- data2 %>% 
#   filter(!(sp == "Taxodium distichum")) %>% # reduced to 4 mort stems
#   filter(!(gr_raw > 5)) %>%  #remove unlikely growth
#   mutate(dbh_raw = exp(dbh)) %>% 
#   filter(dbh_raw < 90)#remove unlikely dbh

# save(data2, file ="R Data/filtered_fit.stan.rdata")
# RUN BELOW -----
load("rdata/filtered_fit.stan.rdata")

# data2 %>% 
#   ggplot(aes(exp(dbh), gr_raw)) +
#   geom_point() +
#   facet_wrap(~sp) +
#   labs(x = "dbh (cm)",
#        y = "growth (cm yr-1)") +
#   mytheme

# plot(density(data$gr^(1/0.47)))
# plot(density(exp(data$dbh)))

# Standardise the explanatory variables:
# **************************************
data_std <- data2
data_std$dbh <- as.numeric(scale(data_std$dbh))
data_std$gr <- as.numeric(scale(data_std$gr))
data = data2

# Add the predictors and the intercept in a matrix format:
compose <- tidybayes::compose_data(data_std)
names(compose)
compose$x <- matrix(c(rep(1, nrow(data_std)), 
                      data_std$dbh,
                      data_std$gr), ncol = 3)
compose$K <- ncol(compose$x)
# 
# MORTALITY PROBABILITY estimates; Careful - this might be very lazy and requires a lot of memory; use the RData!  ------
# Build the survival model - see corresponding stan file --

# Sampler parameters
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
n_cores <- n_chains
seed <- 2024
# For control = list(adapt_delta = ..., max_treedepth = ...):
delta <- 0.95
treedepth <- 12

mod2 <- stan(file = "codes/stan_model_for_mortality_probability.stan", 
            data = compose, 
            chains = n_chains, cores = n_cores, iter = n_iter, 
            warmup = n_warmup, seed = seed, 
            control = list(adapt_delta = delta, max_treedepth = treedepth),
            include = FALSE, pars = c("z_u_sp", "a_plot_z", "a_interval_z"),
            init = 0) # does the trick to avoid getting "rejecting initial values" message problems


# Run this only once, to not have to use "recover_type()" when using spread_draws:
mod2 %<>% recover_types(data) # retrieve original sp and plot names

# load fitted model 
load("rdata/stan_mort_mod_output_nc_mod2_cov_matrix.rdata")
mod = mod2

# Summary of the model output:
# using spread_draws or gather_draws (tidybayes package):
mod %>% 
  tidybayes::spread_draws(a, sigma_a_sp)%>%
  mutate(a = inv.logit(a)) %>% 
  ggplot(aes(x = a)) +
  geom_density() +
  labs(x = "Death probability (yr-1)") +
  mytheme

# hyperparameters and group-level parameters:
mod %>%
  tidybayes::spread_draws(a, sigma_a_sp, a_sp[sp]) %>%
  mutate(a_sp_abs = a + a_sp) %>% 
  mutate(a_sp_abs = inv.logit(a_sp_abs))

post <- mod2 %>%
  tidybayes::spread_draws(a, b_dbh, b_gr,
                          u_sp[sp, param]) %>% # u_sp contains sp-level (row) param. deviations (columns)
  # Transform the column numbers of u_sp (1, 2, and, 3) into what they correspond to:
  rowwise %>% 
  mutate(param = ifelse(param == 1, "a_sp",
                        ifelse(param == 2, "b_dbh_sp", "b_gr_sp"))) %>% 
  ungroup %>% 
  # Pivot to wider format to have each sp-level parameter deviation as a separate column:
  pivot_wider(names_from = param, values_from = u_sp) %>% 
  # Calculate sp-level parameters (summing the hyperparameters to their respective sp-level deviations):
  mutate(spcoef_a = a + a_sp,
         spcoef_b_dbh = b_dbh + b_dbh_sp,
         spcoef_b_gr = b_gr + b_gr_sp) %>%
  # remove the sp-level deviation values (not necessary anymore)
  dplyr::select(-c(a_sp, b_dbh_sp, b_gr_sp))
# save(post, file = "rdata/post.mod.all_cov_matrix.rdata")


#  MORTALITY PROBABILITY estimates: model for early developed forests; Careful - this might be very lazy and requires a lot of memory; use the RData" -----
rm(list = ls())
# load("rdata/data_early_succ_25.RData")
data_std <- data_early
data_std$dbh <- as.numeric(scale(data_std$dbh))
data_std$gr <- as.numeric(scale(data_std$gr))

compose <- tidybayes::compose_data(data_std)

# Sampler parameters
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
n_cores <- n_chains
seed <- 2024
# For control = list(adapt_delta = ..., max_treedepth = ...):
delta <- 0.95
treedepth <- 12

mod_early <- stan(file = "codes/stan_model_for_mortality_probability.stan",
                  data = compose, 
                  chains = n_chains, cores = n_cores, iter = n_iter, 
                  warmup = n_warmup, seed = seed, 
                  control = list(adapt_delta = delta, max_treedepth = treedepth),
                  include = FALSE, pars = c("z_u_sp", "a_plot_z", "a_interval_z"),
                  init = 0) # does the trick to avoid getting "rejecting initial values" message problems

mod_early %<>% recover_types(data_early) # retrieve original sp and plot names
#save(mod_early, file="rdata/stan_mort_mod_output_nc_mod_early_cov_matrix_25.rdata")
# load fitted model 
load("rdata/stan_mort_mod_output_nc_mod_early_cov_matrix_25.rdata")
mod <- mod_early

# Summary of the model output using spread_draws or gather_draws (tidybayes package):
mod %>% 
  tidybayes::spread_draws(a, sigma_a_sp)%>%
  mutate(a = inv.logit(a)) %>% 
  ggplot(aes(x = a)) +
  geom_density() +
  labs(x = "Death probability (yr-1)") +
  mytheme

# hyperparameters and group-level parameters:
mod %>%
  # recover_types(data) %>%
  tidybayes::spread_draws(a, sigma_a_sp, a_sp[sp]) %>%
  mutate(a_sp_abs = a + a_sp) %>% 
  mutate(a_sp_abs = inv.logit(a_sp_abs))

post_early <- mod %>%
  tidybayes::spread_draws(a, b_dbh, b_gr,
                          u_sp[sp, param]) %>% # u_sp contains sp-level (row) param. deviations (columns)
  # Transform the column numbers of u_sp (1, 2, and, 3) into what they correspond to:
  rowwise %>% 
  mutate(param = ifelse(param == 1, "a_sp", ifelse(param == 2, "b_dbh_sp", "b_gr_sp"))) %>% 
  ungroup %>% 
  # Pivot to wider format to have each sp-level parameter deviation as a separate column:
  pivot_wider(names_from = param, values_from = u_sp) %>% 
  # Calculate sp-level parameters (summing the hyperparameters to their respective sp-level deviations):
  mutate(spcoef_a = a + a_sp,
         spcoef_b_dbh = b_dbh + b_dbh_sp,
         spcoef_b_gr = b_gr + b_gr_sp) %>%
  # remove the sp-level deviation values (not necessary anymore)
  dplyr::select(-c(a_sp, b_dbh_sp, b_gr_sp))
#save(post_early, file = "rdata/post.mod.early_cov_matrix_25.rdata")

# MORTALITY PROBABILITY estimates: model for late developed forests; Careful - this might be very lazy and requires a lot of memory; use the RData" ----
rm(list = ls())
#load("rdata/data_late_succ_75.RData")
data_std <- data_late
data_std$dbh <- as.numeric(scale(data_std$dbh))
data_std$gr <- as.numeric(scale(data_std$gr))

compose <- tidybayes::compose_data(data_std)

# Add the predictors and the intercept in a matrix format:
compose$x <- matrix(c(rep(1, nrow(data_std)), 
                      data_std$dbh,
                      data_std$gr), ncol = 3)
compose$K <- ncol(compose$x)

# Sampler parameters
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
n_cores <- n_chains
seed <- 2024
# For control = list(adapt_delta = ..., max_treedepth = ...):
delta <- 0.95
treedepth <- 12

mod_late <- stan(file = "codes/stan_model_for_mortality_probability.stan", 
                 data = compose, 
                 chains = n_chains, cores = n_cores, iter = n_iter, 
                 warmup = n_warmup, seed = seed, 
                 control = list(adapt_delta = delta, max_treedepth = treedepth),
                 include = FALSE, pars = c("z_u_sp", "a_plot_z", "a_interval_z"),
                 init = 0) # does the trick to avoid getting "rejecting initial values" message problems

mod_late %<>% recover_types(data_late) # retrieve original sp and plot names
# save(mod_late, file = "rdata/stan_mort_mod_output_nc_mod_late_cov_matrix_75.rdata")

# load fitted model 
load("rdata/post.mod.late_cov_matrix_75.rdata")
mod <- mod_late

# Summary of the model output using spread_draws or gather_draws (tidybayes package):
mod %>% 
  tidybayes::spread_draws(a, sigma_a_sp)%>%
  mutate(a = inv.logit(a)) %>% 
  ggplot(aes(x = a)) +
  geom_density() +
  labs(x = "Death probability (yr-1)") +
  mytheme

# hyperparameters and group-level parameters:
mod %>%
  # recover_types(data) %>%
  tidybayes::spread_draws(a, sigma_a_sp, a_sp[sp]) %>%
  mutate(a_sp_abs = a + a_sp) %>% 
  mutate(a_sp_abs = inv.logit(a_sp_abs))

post_late <- mod %>%
  tidybayes::spread_draws(a, b_dbh, b_gr,u_sp[sp, param]) %>% # u_sp contains sp-level (row) param. deviations (columns)
  # Transform the column numbers of u_sp (1, 2, and, 3) into what they correspond to:
  rowwise %>% 
  mutate(param = ifelse(param == 1, "a_sp", ifelse(param == 2, "b_dbh_sp", "b_gr_sp"))) %>% 
  ungroup %>% 
  # Pivot to wider format to have each sp-level parameter deviation as a separate column:
  pivot_wider(names_from = param, values_from = u_sp) %>% 
  # Calculate sp-level parameters (summing the hyperparameters to their respective sp-level deviations):
  mutate(spcoef_a = a + a_sp,
         spcoef_b_dbh = b_dbh + b_dbh_sp,
         spcoef_b_gr = b_gr + b_gr_sp) %>%
  # remove the sp-level deviation values (not necessary anymore)
  dplyr::select(-c(a_sp, b_dbh_sp, b_gr_sp))
#save(post_late, file = "rdata/post.mod.late_cov_matrix_75.rdata")

# final rdata with all posteriors ----
load("rdata/post.mod.all_cov_matrix.rdata")
load("rdata/post.mod.late_cov_matrix_75.rdata")
load("rdata/post.mod.early_cov_matrix_25.rdata")
