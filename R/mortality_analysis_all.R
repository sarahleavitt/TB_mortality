#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program performs the mortality Bayesian meta-analysis for all studies
# It takes about 40 minutes to run
##############################################################################

options(scipen=999)
options(digits = 10)
set.seed(150183)

rm(list = ls())
source("R/utils.R")
reload_source()
source("R/mortality_functions.R")

#Reading in individual mortality data
mortality <- read.csv("data/mortality_data.csv")


#### Set-up ----------------------------------------------------------------------------------------

#Creating interval variables
#interval = 1 implies x1 < t <= x2
#dinterval() isn't working as it should so make all obs have interval = 1
#and set x2 to be 10000 (close enough to infinity) for right censored
mortality <- mortality %>% 
  mutate(time = ifelse(death_tb == 0, NA, interval_l),
         interval = 1,
         x1 = interval_l,
         x2 = ifelse(death_tb == 0, 10000, interval_r))

# n.iter <- 31000
# n.burnin <- 1000
# n.thin <- 30

n.iter <- 100
n.burnin <- 10
n.thin <- 1


#Subsetting and formatting data
mortality_strata <- mortality %>%
  #Removing severity stratified mortality data for study 79_1023 (only 4-year follow-up)
  filter(!cohort_id %in% c("79_1023_5", "79_1023_6", "79_1023_7")) %>%
  mutate(study_id_num = as.numeric(factor(study_id)),
         pre1930 = as.numeric(time_period == "pre-1930"),
         northamerica = as.numeric(location == "North America"),
         sanatorium = as.numeric(sanatorium == "Yes"))

#Getting information for each study
study_data <- mortality_strata %>%
  group_by(study_id_num) %>%
  summarize(pre1930 = first(pre1930),
            northamerica = first(northamerica),
            sanatorium = first(sanatorium))


#### Running Models --------------------------------------------------------------------------------

# Complete Model ####

output_comp <- run_comp(mortality_strata, n.iter = n.iter, n.burnin = n.burnin,
                        n.thin = n.thin)


#### Fixed effect for strata separately ####

output_time <- run_one_strata(mortality_strata, study_data, "pre1930",
                              n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)

output_loc <- run_one_strata(mortality_strata, study_data, "northamerica",
                              n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)

output_san <- run_one_strata(mortality_strata, study_data, "sanatorium",
                              n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)


#### Fixed effect for strata combined ####

output_all <- run_all_strata(mortality_strata, study_data, n.iter = n.iter,
                             n.burnin = n.burnin, n.thin = n.thin)



#### Saving Results------------------------------------------------------------------------------

data_comp <- getData(mortality_comp)
res_comp <- output_comp$res
eval_comp <- output_comp$eval


save(res_comp_all, res_sev_all, eval_comp_all, eval_sev_all, data_comp_all, data_sev_all,
     file = "R/bayesian_all.RData")


