#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

#################################################################################
# This program performs the mortality Bayesian meta-analysis for the US studies
# It takes about 30 minutes to run
#################################################################################

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


#### US Studies: Complete Model ####

#Subsetting and formatting data
mortality_comp_us <- mortality %>%
  filter(study_id %in% c("1029", "93", "45", "63", "67", "90_1016")) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)))

#Running model
output_comp_us <- run_comp(mortality_comp_us,
                            n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)



#### US Studies: Stratified Model ####

#Subsetting and formatting data
mortality_sev_us <- mortality %>%
  filter(study_id %in% c("1029", "93", "45", "63", "67", "90_1016"),
         severity != "Unknown") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)),
         sev_mod = as.numeric(severity == "Moderate"),
         sev_adv = as.numeric(severity == "Advanced"),
         sev_unk = as.numeric(severity == "Unknown"))

#Running model
output_sev_us <- run_sev(mortality_sev_us,
                          n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)





#### Saving Results------------------------------------------------------------------------------

data_comp_us <- getData(mortality_comp_us)
data_sev_us <- getData(mortality_sev_us)
res_comp_us <- output_comp_us$res
res_sev_us <- output_sev_us$res
eval_comp_us <- output_comp_us$eval
eval_sev_us <- output_sev_us$eval

save(res_comp_us, res_sev_us, eval_comp_us, eval_sev_us, data_comp_us, data_sev_us,
     file = "R/bayesian_us.RData")


