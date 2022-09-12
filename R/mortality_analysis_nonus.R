#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

####################################################################################
# This program performs the mortality Bayesian meta-analysis for the Non-US studies
# It takes about 25 to run
####################################################################################

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

n.iter <- 31000
n.burnin <- 1000
n.thin <- 30

# n.iter <- 100
# n.burnin <- 10
# n.thin <- 1



#### Running Models --------------------------------------------------------------------------------


#### Non-US Studies: Complete Model ####

#Subsetting and formatting data
mortality_comp_nonus <- mortality %>%
  filter(!study_id %in% c("1029", "93", "45", "63", "67", "90_1016"),
         !cohort_id %in% c("79_1023_5", "79_1023_6", "79_1023_7")) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)))

#Running model
output_comp_nonus <- run_comp(mortality_comp_nonus,
                            n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)



#### Non-US Studies: Stratified Model ####

#Subsetting and formatting data
mortality_sev_nonus <- mortality %>%
  filter(!study_id %in% c("1029", "93", "45", "63", "67", "90_1016"),
         severity != "Unknown") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)),
         sev_mod = as.numeric(severity == "Moderate"),
         sev_adv = as.numeric(severity == "Advanced"),
         sev_unk = as.numeric(severity == "Unknown"))

#Running model
output_sev_nonus <- run_sev(mortality_sev_nonus,
                          n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)





#### Saving Results------------------------------------------------------------------------------

data_comp_nonus <- getData(mortality_comp_nonus)
data_sev_nonus <- getData(mortality_sev_nonus)
res_comp_nonus <- output_comp_nonus$res
res_sev_nonus <- output_sev_nonus$res
eval_comp_nonus <- output_comp_nonus$eval
eval_sev_nonus <- output_sev_nonus$eval

save(res_comp_nonus, res_sev_nonus, eval_comp_nonus, eval_sev_nonus, data_comp_nonus, data_sev_nonus,
     file = "R/bayesian_nonus.RData")


