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

#Reading in individual mortality data and analysis results
mortality <- read.csv("data/mortality_data.csv")
load('R/bayesian_mortality.RData')
load('R/bayesian_all.RData')
load('R/bayesian_us.RData')
load('R/bayesian_nonus.RData')
load('R/bayesian_post.RData')
load('R/bayesian_pre.RData')
load('R/bayesian_san.RData')



#### Table of Main Results -------------------------------------------------------------------------

#Combining raw tables from results lists
raw_tab <- bind_rows(form_comp_all$param, form_sev_all$param,
                     form_comp_us$param, form_sev_us$param,
                     form_comp_nonus$param, form_sev_nonus$param,
                     form_comp_post$param, form_sev_post$param,
                     form_comp_pre$param, form_sev_pre$param,
                     form_san$param, form_nosan$param)

#Adding formatted, ordered labels
raw_tab <- raw_tab %>%
  mutate(Severity = ifelse(is.na(severity), label,
                           ifelse(label == "Severity_all", paste0(severity, "_all"),
                                  ifelse(label == "Severity_us", paste0(severity, "_us"),
                                         ifelse(label == "Severity_nonus", paste0(severity, "_nonus"),
                                                ifelse(label == "Severity_post", paste0(severity, "_post"),
                                                       ifelse(label == "Severity_pre", paste0(severity, "_pre"),
                                                              label)))))),
         Severity = factor(Severity, level = c("Minimal_all", "Moderate_all", "Advanced_all", "Combined_all",
                                               "Minimal_us", "Moderate_us", "Advanced_us", "Combined_us",
                                               "Minimal_nonus", "Moderate_nonus", "Advanced_nonus", "Combined_nonus",
                                               "Minimal_post", "Moderate_post", "Advanced_post", "Combined_post",
                                               "Minimal_pre", "Moderate_pre", "Advanced_pre", "Combined_pre",
                                               "Sanatorium/hospital", "Non-Sanatorium"))) %>%
  arrange(Severity)

#Extracting the survival probabilities
pred1_tab <- raw_tab %>%
  filter(value == "pred1") %>%
  mutate(`1-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                             round(cilb, 2), ", ",
                                             round(ciub, 2), ")")) %>%
  select(Severity, `1-Year Survival (95% CI)`)

pred5_tab <- raw_tab %>%
  filter(value == "pred5") %>%
  mutate(`5-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                             round(cilb, 2), ", ",
                                             round(ciub, 2), ")")) %>%
  select(Severity, `5-Year Survival (95% CI)`)

pred10_tab <- raw_tab %>%
  filter(value == "pred10") %>%
  mutate(`10-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                              round(cilb, 2), ", ",
                                              round(ciub, 2), ")")) %>%
  select(Severity, `10-Year Survival (95% CI)`)


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


#UPDATE

#Finding number of papers, cohorts, individual for each analysis
data_comp_all2 <- as.data.frame(t(data_comp_all[[2]]))
data_comp_all2$Severity <- "Combined_all"
data_comp_us2 <- as.data.frame(t(data_comp_us[[2]]))
data_comp_us2$Severity <- "Combined_us"
data_comp_nonus2 <- as.data.frame(t(data_comp_nonus[[2]]))
data_comp_nonus2$Severity <- "Combined_nonus"
data_comp_post2 <- as.data.frame(t(data_comp_post[[2]]))
data_comp_post2$Severity <- "Combined_post"
data_comp_pre2 <- as.data.frame(t(data_comp_pre[[2]]))
data_comp_pre2$Severity <- "Combined_pre"
data_san2 <- as.data.frame(t(data_san[[2]]))
data_san2$Severity <- "Sanatorium/hospital"
data_nosan2 <- as.data.frame(t(data_nosan[[2]]))
data_nosan2$Severity <- "Non-Sanatorium"

counts_comp <- bind_rows(data_comp_all2, data_comp_us2, data_comp_nonus2,
                         data_comp_post2, data_comp_pre2, data_san2, data_nosan2)

#Counts stratified by severity for all studies
mortality_sev <- mortality %>% filter(severity != "Unknown")

counts_all <- mortality_sev %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  mutate(Severity = paste0(severity, "_all"))


#Counts stratified by severity for US studies
mortality_us <- mortality %>% filter(severity != "Unknown",
                                     study_id %in% c("1029", "93", "45", "63", "67", "90_1016"))

counts_us <- mortality_us %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  mutate(Severity = paste0(severity, "_us"))


#Counts stratified by severity for Non-US studies
mortality_nonus <- mortality %>% filter(severity != "Unknown",
                                     !study_id %in% c("1029", "93", "45", "63", "67", "90_1016"))

counts_nonus <- mortality_nonus %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  mutate(Severity = paste0(severity, "_nonus"))


#Counts stratified by severity for US post-1930s subset
mortality_post <- mortality %>% filter(severity != "Unknown", study_id %in% c("1029", "93", "45"))

counts_post <- mortality_post %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  mutate(Severity = paste0(severity, "_post"))


#Counts stratified by severity for US pre-1930s subset
mortality_pre <- mortality %>% filter(severity != "Unknown", study_id %in% c("63", "67", "90_1016"))

counts_pre <- mortality_pre %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  mutate(Severity = paste0(severity, "_pre"))


counts_initial <- bind_rows(counts_comp, counts_all, counts_us, counts_nonus,
                            counts_post, counts_pre) %>%
  select(Severity, `Number of Studies` = nStudies,
         `Number of Cohorts` = nCohorts,
         `Number of Individuals` = nIndividuals)



#Combining the tables
final_tab <- dist_tab %>%
  full_join(pred1_tab, by = c("Severity")) %>%
  full_join(pred5_tab, by = c("Severity")) %>%
  full_join(pred10_tab, by = c("Severity")) %>%
  full_join(counts_initial, by = c("Severity"))


#Separating tables for paper
all_tab <- final_tab %>% filter(grepl("_all", Severity))
us_tab <- final_tab %>% filter(grepl("_us", Severity))
nonus_tab <- final_tab %>% filter(grepl("_nonus", Severity))
post_tab <- final_tab %>% filter(grepl("_post", Severity))
pre_tab <- final_tab %>% filter(grepl("_pre", Severity))
san_tab <- final_tab %>% filter(grepl("San", Severity))

#Variance of frailty terms
theta <- raw_tab %>%
  filter(value == "theta") %>%
  mutate(theta = round(est, 2)) %>%
  select(label, theta)

