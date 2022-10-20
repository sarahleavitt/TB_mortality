#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program creates figures and tables of the mortality analysis results
##############################################################################

options(scipen=999)
options(digits = 10)

rm(list = ls())
source("R/utils.R")
reload_source()

#Reading in the study_id correspondence table
studyid <- read.csv("data/study_id.csv")

load('R/bayesian_raw.RData')
load('R/bayesian_clean.RData')



#### Table of Main Results -----------------------------------------------------

#Combining raw tables from results lists
raw_tab <- bind_rows(form_comb$param, form_pre$param, form_post$param,
                     form_namerica$param, form_europe$param, form_yessan$param,
                     form_nosan$param)

#Extracting the survival probabilities
pred1_tab <- raw_tab %>%
  filter(value == "pred1") %>%
  mutate(`1-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                             round(cilb, 2), ", ",
                                             round(ciub, 2), ")")) %>%
  select(label, `1-Year Survival (95% CI)`)

pred5_tab <- raw_tab %>%
  filter(value == "pred5") %>%
  mutate(`5-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                             round(cilb, 2), ", ",
                                             round(ciub, 2), ")")) %>%
  select(label, `5-Year Survival (95% CI)`)

pred10_tab <- raw_tab %>%
  filter(value == "pred10") %>%
  mutate(`10-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                              round(cilb, 2), ", ",
                                              round(ciub, 2), ")")) %>%
  select(label, `10-Year Survival (95% CI)`)


#Extracting median survival
med_tab <- raw_tab %>%
  filter(value == "median") %>%
  mutate(`Median Survival (95% CI)` = paste0(round(est, 2), " (",
                                              round(cilb, 2), ", ",
                                              round(ciub, 2), ")")) %>%
  select(label, `Median Survival (95% CI)`)


#Extracting the distribution parameters
sdlog <- raw_tab %>%
  filter(value == "sdlog") %>%
  select(label, sdlog = est)

dist_tab <- raw_tab %>%
  filter(value == "meanlog") %>%
  full_join(sdlog, by = "label") %>%
  mutate(`Survival Distribution` = paste0("lognormal(", round(est, 2), ", ",
                                          round(sdlog, 2), ")")) %>%
  select(label, `Survival Distribution`)


# Finding number of papers, cohorts, individuals
data_comb2 <- as.data.frame(t(data_comb[[2]]))
data_comb2$label <- "Combined"
data_pre2 <- as.data.frame(t(data_pre[[2]]))
data_pre2$label <- "Pre-1930"
data_post2 <- as.data.frame(t(data_post[[2]]))
data_post2$label <- "Post-1930"
data_namerica2 <- as.data.frame(t(data_namerica[[2]]))
data_namerica2$label <- "North America"
data_europe2 <- as.data.frame(t(data_europe[[2]]))
data_europe2$label <- "Europe"
data_yessan2 <- as.data.frame(t(data_yessan[[2]]))
data_yessan2$label <- "Sanatorium/hospital"
data_nosan2 <- as.data.frame(t(data_nosan[[2]]))
data_nosan2$label <- "Non-Sanatorium"

counts <- bind_rows(data_comb2, data_pre2, data_post2, data_namerica2,
                    data_europe2, data_yessan2, data_nosan2) %>%
  select(label, `Number of Studies` = nStudies,
         `Number of Cohorts` = nCohorts,
         `Number of Individuals` = nIndividuals)


#Combining the tables
final_tab <- dist_tab %>%
  full_join(pred1_tab, by = c("label")) %>%
  full_join(pred5_tab, by = c("label")) %>%
  full_join(pred10_tab, by = c("label")) %>%
  full_join(med_tab, by = c("label")) %>%
  full_join(counts, by = c("label"))

write.csv(final_tab, "data/mortality_results_table.csv")


#Variance of frailty terms
theta <- raw_tab %>%
  filter(value == "theta") %>%
  mutate(theta = round(est, 2)) %>%
  select(label, theta)

