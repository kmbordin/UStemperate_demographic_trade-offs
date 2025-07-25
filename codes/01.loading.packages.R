########################################################################
####### PACKAGES AND FUNCTIONS #########################################
########################################################################
rm(list = ls())

# loading packages -----
library(tidyverse)  # tidy data
library(data.table) # for data organisation
library(reshape2)   # for data organisation 
library(here)       # for loading data
library(doBy)       # to summarise data
library(smatr)		  # for MA/SMA regression
library(boot)       # for regressions
library(ggpubr)     # for plots
library(tidyr)      # unnest lists
library(patchwork)  # plots organisation
library(ggpmisc)    # SMA model fit
library(rstan)      # stan code to estimate mortality
library(tidybayes)  # bayesian results in tidy format
library(ggridges)   # bayesian plots
library(magrittr)   # pipes and tidy
library(modelr)     # extract bayesian results
library(ggforce)    # figures
library(dplyr)      # data manipulation
set.seed(2024)      # turn on to generate the same result in permutations


# theme for plots ----
mytheme <- theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 17),
        plot.caption = element_text(size = 13))
my_theme <- theme_light() +
  theme (panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position = "none",
         aspect.ratio = 1,
         text = element_text(size=10),
         plot.margin=unit(c(0, 0, 0 , 0), "cm"))


# filter species to run mortality rates ----
species_list_to_mortality <- function (x, min.n.surv, min.n.mort) {
  # x = data
  # min.n.surv = set the minimum number of survivors
  # min.n.mort = set the minimum number of deads
  
  # info of surviving stems
  splist.surv <- x %>%
    filter(DFstatus3 ==0) %>% # stem status; zero is alive
    count(sp) %>%
    filter(n >= min.n.surv) # filter only the species that meet the threshold
  
  # info of dead stems
  splist.mort <- x %>%
    filter(DFstatus3 ==1) %>% # stem status; one is dead
    count(sp) %>%
    filter(n >= min.n.mort) # filter only the species that meet the threshold
  
  # species list from live stems
  splist <- splist.surv %>%
    filter(sp %in% splist.mort$sp) # filter the same species across dead and alive trees
  
  # filtering the species information from the dataset
  x <- x %>%
    filter(sp %in% splist$sp)
  
  # creating a new table with the filtered species
  table <- data.frame(tag = x$tag, # tag info
                      sp = x$sp, # species info
                      plot = x$plot.id, # plot id info
                      mort = x$DFstatus3,# stem status; maintain the code of 1 for dead stems and 0 for live
                      dbh2 = x$dbh2,# previous dbh (census 2)
                      time.growth = x$date2-x$date1) %>% # time interval from census 1 to 2
    mutate(prev.growth = (x$dbh2-x$dbh1)/time.growth, # previous growth (census2-census1)
           time = x$date3-x$date2,#time census interval between census 2 and 3
           gr = ifelse(test= prev.growth>=0, yes= prev.growth^0.47, no= prev.growth), # to reduce skewness
           dbh = log(dbh2)) %>%
    filter(sp != "Indet indet") # remove unidentified species
  
  table$gr[table$gr < 0] = -(abs(table$gr[table$gr < 0])^0.47) #reduce skewness of non-positive growth
  
  table <- table %>% 
    select(tag, sp, mort, dbh, gr, time, plot) %>% #columns to be modeled in stan to obtain the mortality probabilities 
    mutate(tag = as.factor(tag), #ensuring the right structure of data
           sp = as.character(sp),
           plot = as.factor(plot))
}
# bootstrap growth rates ----
obtain_sd_growth95 <- function(x, nboot = 999){
  # x = dataset containing tree growth rates
  # growth.dbh = the column with the real growth
  growth95.boot <- c() #save results
  obs.growth95 <- quantile(x$growth.dbh, probs = c(0.95), na.rm=TRUE) # 95th quantile of growth
  
  for(i in 1:nboot){
    temp <- sample(1:nrow(x), nrow(x), replace = T) # sampling with replace
    temp2 <- x[temp, "growth.dbh"] 
    growth95.sample.i <- quantile(temp2, probs = c(0.95), na.rm=TRUE) # sample the 95th quantiles
    growth95.boot <- c(growth95.boot, growth95.sample.i) # save the results in a vector
  }
  c(obs.growth95 = growth95.sample.i, # return the values
    sd = sd(growth95.boot), 
    quantile(growth95.boot, probs = c(0.025, 0.975)))
}

# detect skewness across the data - this function in available in the CTFSRPackage ----- 
# http://ctfs.si.edu/ctfsdev/CTFSRPackageNew/
skewness <- function (x) {
  x = x[!is.na(x)]
  n = length(x)
  if (length(x) <= 2) 
    return(0)
  mn = mean(x)
  sumcube = sum((x - mn)^3)
  sumsq = sum((x - mn)^2)
  biased = sqrt(n) * sumcube/(sumsq^1.5)
  correction = sqrt(n) * sqrt((n - 1)/(n - 2))
  return(biased * correction)
}
# extract model posteriors ----
preds <-  function (data, post){
  # data = dataset used to generate the estimates
  # post = posteriors from the models
  data_std = data
  data_std$dbh <- as.numeric(scale(data_std$dbh))
  data_std$gr <- as.numeric(scale(data_std$gr))
  
  # Build a grid of covariate value combinations for which we want posterior predictions:
  # adequate accordingly
  dbh_val <- (log(12.7) - mean(data$dbh)) / sd(data$dbh) # corresponds to dbh of 12.7 cm
  gr_val <- (0^(0.47) - mean(data$gr)) / sd(data$gr) # corresponds to gr of 0 cm yr-1
  
  grid <- modelr::data_grid(data = data_std,
                            gr = gr_val,
                            dbh = dbh_val,
                            sp)
  
  # Predictions from the fitted model only:
  predict <- post %>% 
    # add grid of covariate values to the posterior draws:
    left_join(grid) %>% 
    # Prediction per species (using species-level parameters (hyperparam + sp deviation) only):
    mutate(pred_sp = spcoef_a + spcoef_b_dbh * dbh + spcoef_b_gr * gr) %>% 
    rename(pred_sp_m = pred_sp) %>% 
    mutate(pred_sp_m_invl = inv.logit(pred_sp_m)) %>% # inverse logit to probability scale
    # back-transform covariate, if they were standardised prior to running the model:
    mutate(dbh_raw = (dbh * sd(data$dbh)) + mean(data$dbh),
           dbh_raw = exp(dbh_raw),
           gr_raw = (gr * sd(data$gr)) + mean(data$gr),
           gr_raw = gr_raw^(1/0.47))
  # Summarise per species:
  (predict_summ <- predict %>% 
      group_by(sp, dbh_raw, gr_raw, dbh, gr) %>% 
      point_interval(pred_sp_m_invl, 
                     .point = mean, .interval = qi, .width = 0.9) %>% 
      ungroup)
  
} 
# plotting model posteriors ----
plot_post <- function(pred,var){
  # pred = predictions extracted from the model posteriors
  # var = corresponding variable
  pred %>% 
    ggplot(aes(pred_sp_m_invl, sp)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper)) +
    mytheme +
    scale_y_discrete(limits=rev)+
    labs(x = "Mortality probability (yr-¹)",
         subtitle = "Species level mortality probability\n at zero growth and 12.7 cm of diameter",
         title = var) +
    # coord_cartesian(xlim = c(0, 0.15)) + # to limit the x-axis range
    theme(axis.text.y = element_text(size = 12, face = "italic"),
          axis.title.y = element_blank(),
          plot.subtitle = element_text(hjust = 0.5),
          plot.title = element_text(hjust = 0.5))
}