#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program compares the results from the different analyses for the 
# mortality paper
##############################################################################

options(scipen=999)
options(digits = 10)

#rm(list = ls())
source("R/utils.R")
reload_source()

#Reading in the study_id correspondence table
studyid <- read.csv("data/study_id.csv")

#Reading in individual mortality data and analysis results
mortality <- read.csv("data/mortality_data.csv")
load('R/bayesian_comp.RData')
load('R/bayesian_onestrata.RData')
load('R/bayesian_allstrata.RData')
load('R/bayesian_separate.RData')


#Combining raw tables from results lists
raw_tab <- bind_rows(form_comp$param, form_all,
                     form_time$param, form_loc$param, form_san$param,
                     form_pre$param, form_post$param, 
                     form_namerica$param, form_europe$param,
                     form_yessan$param, form_nosan$param)

#Extracting the survival probabilities
pred1_tab <- raw_tab %>%
  filter(value == "pred1") %>%
  mutate(`1-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                             round(cilb, 2), ", ",
                                             round(ciub, 2), ")")) %>%
  select(label, strata, `1-Year Survival (95% CI)`)

pred10_tab <- raw_tab %>%
  filter(value == "pred10") %>%
  mutate(`10-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                              round(cilb, 2), ", ",
                                              round(ciub, 2), ")")) %>%
  select(label, strata, `10-Year Survival (95% CI)`)

med_tab <- raw_tab %>%
  filter(value == "median") %>%
  mutate(`Median Survival Time (95% CI)` = paste0(round(est, 2), " (",
                                                  round(cilb, 2), ", ",
                                                  round(ciub, 2), ")")) %>%
  select(label, strata, `Median Survival Time (95% CI)`)


#Extracting the distribution parameters
sdlog <- raw_tab %>%
  filter(value == "sdlog") %>%
  select(label, sdlog = est)

dist_tab <- raw_tab %>%
  filter(value == "meanlog") %>%
  full_join(sdlog, by = "label") %>%
  mutate(`Survival Distribution` = paste0("lognormal(", round(est, 2), ", ",
                                          round(sdlog, 2), ")")) %>%
  select(Severity, `Survival Distribution`)