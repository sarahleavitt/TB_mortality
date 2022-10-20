#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program transforms the life-table study data extracted from the 
# publications and transforms it into individual-level data
##############################################################################

options(scipen=999)
options(digits = 10)

rm(list = ls())
source("R/utils.R")
reload_source()

## Reading in all of the datasets
dataList <- read_excel_allsheets("data/pre_chemo_data.xlsx")
dataList$`Data dictionary` <- NULL

## Reading in study information
studyid <- read.csv("data/study_id.csv") %>%
  select(-start_type)


#### Overall counts -------------------------------------------------------------------------------

#Removing the severity data for 75_23 because it is the same people as the full study data in 75_1023
countList <- dataList[!names(dataList) %in% c("79_1023_sev")]

pull_first_row <- function(paper){
  
  first_row <- paper %>%
    group_by(cohort_id) %>%
    arrange(interval_l) %>%
    select(study_id, paper_id, cohort_id, n, c1a, c2) %>%
    mutate(paper_id = as.character(paper_id),
           study_id = as.character(study_id)) %>%
    slice(1)
  
  return(first_row)
}

cohorts <- map_dfr(countList, pull_first_row)

#Number of papers
length(unique(cohorts$paper_id))

#Number of studies
length(unique(cohorts$study_id))

#Number of cohorts
length(unique(cohorts$cohort_id))

#Number of patients
sum(cohorts$n)

#Subset to TB mortality data
mortality <- cohorts %>% filter(!is.na(c1a)) 
length(unique(mortality$study_id))
length(unique(mortality$cohort_id))

#Subset to cure data
cure <- cohorts %>%
  filter(!is.na(c2) | study_id == "5_1047") %>%
  bind_rows(select(dataList$`79_1023_sev`, study_id, paper_id, cohort_id, n, c1a, c2))
length(unique(cure$study_id))
length(unique(cure$cohort_id))



#### Mortality Data -------------------------------------------------------------------------------

mortalityList <- dataList[!names(dataList) %in% c("65", "91", "75_1019", "75_1020", "75_1021",
                                                  "75_1019", "75_1019", "75_1021", "79_1054")]

#Converting the study data to individual mortality data
indAll <- map_dfr(mortalityList, studyToInd, outcome = "mortality")

#Formatting individual dataset
mortalityData <- indAll %>%
  mutate(#People who are censored at the study start
         no_survival = death == 0 & interval_l == 0,
         start_type = ifelse(start_type == "Unknown", "Entry", start_type),
         #Adding overall average length of stay to bring start_type = 'Exit' start_type = 'Entry'
         interval_l = ifelse(start_type == "Exit", interval_l + 165/365, interval_l),
         interval_r = ifelse(start_type == "Exit", interval_r + 165/365, interval_r),
         severity = ifelse(severity %in% c("Incipient disease",
                                                 "Minimal disease",
                                                 "Stage I",
                                                 "TB plus I (early)"), "Minimal",
                           ifelse(severity %in% c("Moderately advanced disease",
                                                        "Stage II",
                                                        "TB plus II (intermediate)"), "Moderate",
                                  ifelse(severity %in% c("Far advanced disease",
                                                               "Stage III",
                                                               "TB plus III (advanced)"), "Advanced",
                                         "Unknown"))),
         severity = factor(severity, levels = c("Minimal", "Moderate", "Advanced", "Unknown"))) %>%
  unite(study_sev, study_id, severity, remove = FALSE) %>%
  left_join(studyid, by = "study_id") %>%
  arrange(study_sev)

write.csv(mortalityData, "data/mortality_data.csv", row.names = FALSE)

#Baseline sample size for mortality analysis
nrow(mortalityData)



#### Cure Data with severity ----------------------------------------------------------------------

cureList <- list(dataList$`1029_1055`, dataList$`1029_1056`, dataList$`48_1000_1029`,
                 dataList$`45`, dataList$`67`)

## Converting the study data to individual mortality data
cureData1 <- map_dfr(cureList, studyToInd, outcome = "cure", timepoints = 3)

## Extracting four year survival from 79_1023
cureData2 <- studyToInd(dataList$`79_1023_sev`, outcome = "cure", timepoints = 4)

#Combining Info
cureData <- bind_rows(cureData1, cureData2) %>%
  left_join(studyid, by = "study_id")

write.csv(cureData, "data/cure_data.csv", row.names = FALSE)




#### All Cure Data ----------------------------------------------------------------------------

cureListAll <- list(dataList$`75_1019`,dataList$`75_1020`,dataList$`75_1021`,dataList$`79_1054`,
                    dataList$`5_1047`,dataList$`94`, dataList$`65`,dataList$`91`,
                    dataList$`1029_1055`,dataList$`1029_1056`,dataList$`48_1000_1029`,
                    dataList$`79_1023_sev`,dataList$`45`,dataList$`67`)

cureListAll2 <- lapply(cureListAll,dblTochr)

cureDataAll <- map_dfr(cureListAll2,calcCureRate)

write.csv(cureDataAll, "data/cure_data_all.csv", row.names = FALSE)

