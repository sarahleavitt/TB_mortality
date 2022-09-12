#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program performs the mortality Bayesian meta-analysis for all studies
# It takes about 15 minutes to run
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


#Subsetting and formatting data
mortality_comp_all <- mortality %>%
  #Removing severity stratified mortality data for study 79_1023 (only 4-year follow-up)
  filter(!cohort_id %in% c("79_1023_5", "79_1023_6", "79_1023_7")) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)))

n.iter <- 31000
n.burnin <- 1000
n.thin <- 30

# n.iter <- 100
# n.burnin <- 10
# n.thin <- 1



#### Running Models --------------------------------------------------------------------------------


#### Sanatorium/hospital ####

#Subsetting data
san <- mortality_comp_all %>%
  filter(sanatorium == "Yes") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)))

#Running model
output_san <- run_comp(san)


#### Non-sanatorium ####

#Subsetting data
nosan <- mortality_comp_all %>%
  filter(sanatorium == "No") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)))

#Running model
output_nosan <- run_comp(nosan)




#### Saving Results------------------------------------------------------------------------------

data_san <- getData(san)
data_nosan <- getData(nosan)
res_san <- output_san$res
res_nosan <- output_nosan$res
eval_san <- output_san$eval
eval_nosan <- output_nosan$eval

save(res_san, res_nosan, eval_san, eval_nosan, data_san, data_nosan,
     file = "R/bayesian_san.RData")


