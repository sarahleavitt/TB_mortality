
## Function to load packages
reload_source <- function(){
  if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
  if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
  if (!require('purrr')) install.packages('purrr'); library('purrr')
  if (!require('stringr')) install.packages('stringr'); library('stringr')
  if (!require('naniar')) install.packages('naniar'); library('naniar')
  if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
  if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
  if (!require('readxl')) install.packages('readxl'); library('readxl')
  if (!require('R2jags')) install.packages('R2jags'); library('R2jags')
  if (!require('lattice')) install.packages('lattice'); library('lattice')
  if (!require('mcmcplots')) install.packages('mcmcplots'); library('mcmcplots')
}


## Function to read all sheets from an excel file
read_excel_allsheets <- function(filename, tibble = TRUE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


## Function to turn one sheet of study data into individual-level data
studyToInd <- function(DF, outcome = c("mortality", "cure"), timepoints = c(3, 5, 10)){
  
  #Creating individual-level mortality data
  if(outcome == "mortality"){
    
    
    #Finding all-cause mortality
    if(is.na(DF$c1a_plus_b)[1]){
      DF <- DF %>% mutate(c1a_plus_b = c1a + c1b) 
    }
    
    cohorts <- DF %>%
      group_by(cohort_id) %>%
      group_modify(~ cohortToInd(.x)) %>%
      ungroup()
  }

  #Creating individual-level cure data
  if(outcome == "cure"){
    cohorts <- DF %>%
      group_by(cohort_id) %>%
      group_modify(~ cohortToIndCure(.x, timepoints = timepoints)) %>%
      ungroup()
  }
  
  cohorts <- as.data.frame(cohorts)
}


## Function to turn one cohort into individual mortality data (run by studyToInd)
cohortToInd <- function(DFc){
  
  #Creating individual data
  #death = 1 means that the subject died and death = 0 means they were censored
  #death_tb = 1 means that the subject died of TB
  #death_other = 1 means the subject died of other causes
  
  ind <- NULL
  for(i in 1:nrow(DFc)){
    row <- DFc[i, ]
    interval_l <- row$interval_l
    interval_r <- row$interval_r
    
    #If there is TB verses all-cause mortality, split the two, if not combine
    if(!is.na(DFc$c1a[1])){
      
      #Death from TB
      death_tb <- cbind.data.frame(rep(interval_l, row$c1a),
                                   rep(interval_r, row$c1a),
                                   rep(1, row$c1a))
      names(death_tb) <- c("interval_l", "interval_r", "death_tb")
      
      #Death from other causes
      death_other <- cbind.data.frame(rep(interval_l, row$c1b),
                                    rep(interval_r, row$c1b),
                                    rep(1, row$c1b))
      names(death_other) <- c("interval_l", "interval_r", "death_other")
      
      #Combining death from both causes
      death <- bind_rows(death_tb, death_other) %>%
        replace_na(list(death_tb = 0, death_other = 0)) %>%
        mutate(death = death_other + death_tb)
      
    }else{
      death <- cbind.data.frame(rep(interval_l, row$c1a_plus_b),
                                rep(interval_r, row$c1a_plus_b),
                                rep(1, row$c1a_plus_b))
      names(death) <- c("interval_l", "interval_r", "death")
    }
    
    #Finding those who were lost
    censor <- cbind.data.frame(rep(interval_l, row$l),
                               rep(interval_r, row$l),
                               rep(0, row$l),
                               rep(0, row$l),
                               rep(0, row$l))
    names(censor) <- c("interval_l", "interval_r", "death_tb", "death_other", "death")
    
    indTemp <- bind_rows(death, censor)
    ind <- bind_rows(ind, indTemp)
  }
  
  #Adding individuals who are censored at the end of the follow-up
  total <- DFc$n[1]
  haveOutcome <- nrow(ind)
  endTime <- DFc$interval_r[nrow(DFc)]
  censor <- cbind.data.frame(rep(endTime, total - haveOutcome), 
                             rep(0, total - haveOutcome),
                             rep(0, total - haveOutcome),
                             rep(0, total - haveOutcome))
  names(censor) <- c("interval_l", "death_tb", "death_other", "death")
  
  #Adding study-level characteristics
  ind <- bind_rows(ind, censor) %>%
    mutate(start_type = DFc$start_type[1],
           severity = DFc$severity[1],
           paper_id = as.character(DFc$paper_id[1]),
           study_id = as.character(DFc$study_id[1]))
  
  return(ind)
}


## Function to turn one cohort into individual cure data (run by studyToInd)
cohortToIndCure <- function(DFc, timepoints = c(3, 4, 5, 10)){
  
  nTotal <- as.numeric(DFc[1, "n"])
  
  ind <- NULL
  for(t in timepoints){
    nCure <- as.numeric(DFc[DFc$interval_r == t, "c2"])
    nLost <- sum(DFc[DFc$interval_r <= t, "l"])
    n <- nTotal - nLost
    
    cure <- cbind.data.frame("time" = t, "cure" = rep(1, nCure))
    censor <- cbind.data.frame("time" = t, "cure" = rep(0, n - nCure))
    ind <- bind_rows(ind, cure, censor)
  }
  
  ind <- ind %>% mutate(start_type = DFc$start_type[1],
                        severity = DFc$severity[1],
                        paper_id = as.character(DFc$paper_id[1]),
                        study_id = as.character(DFc$study_id[1]))
  
  return(ind)
}


## Function to get information about the dataset for each run
getData <- function(data){
  
  #Finding concordance to frailty IDs
  tab <- data %>%
    select(study_id, study_id_num) %>%
    filter(!duplicated(study_id))
  
  #Finding n studies, severity groups, cohorts, individuals
  counts <- c("nStudies" = length(unique(data$study_id_num)),
              "nCohorts" = length(unique(data$cohort_id)),
              "nIndividuals" = nrow(data))
  
  return(list(tab, counts))
}


## Function to format the parameter table from the Bayesian analysis
format_param <- function(res, label, fixed){
  
  #Parameter values
  if(fixed == FALSE){
    
    param <- res[row.names(res) %in% c("mu", "sdlog", "theta", "med_comp",
                                       "pred_comp[1]", "pred_comp[5]", "pred_comp[10]"), ]
    names(param) <- c("cilb", "lowerquant", "est", "upperquant", "ciub")
    param <- param %>% mutate(value = row.names(param),
                              value = ifelse(value == "mu", "meanlog",
                                             ifelse(value == "med_comp", "median",
                                                    gsub("_comp\\[|\\]", "", value))),
                              label = label)
  }else{
    param <- res[row.names(res) %in% c("meanlog_min", "meanlog_mod", "meanlog_adv",
                                       "sdlog", "theta", "med_min", "med_mod", "med_adv",
                                       "pred_min[1]", "pred_min[5]", "pred_min[10]",
                                       "pred_mod[1]", "pred_mod[5]", "pred_mod[10]",
                                       "pred_adv[1]", "pred_adv[5]", "pred_adv[10]"), ]
    names(param) <- c("cilb", "lowerquant", "est", "upperquant", "ciub")
    param <- param %>% mutate(label = label,
                              severity = ifelse(grepl("min", row.names(.)), "Minimal",
                                                ifelse(grepl("mod", row.names(.)), "Moderate",
                                                       ifelse(grepl("adv", row.names(.)), "Advanced", NA))),
                              value = gsub("_[a-z]*$|_[a-z]*\\[|\\]", "", row.names(.)),
                              value = ifelse(value == "med", "median", value))
  }
  return(param)
}

## Function to format the overall survival and density table from the Bayesian analysis
format_surv_dens <- function(param, res, label, fixed){
  
  if(fixed == FALSE){
    
    #Overall survival and density curves
    sdlog <- param %>% filter(value == "sdlog") %>% pull(est)
    meanlog <- param %>% filter(value == "meanlog") %>% pull(est)
    x <- seq(0, 30, 0.01)
    dens <- dlnorm(x, meanlog, sdlog)
    surv <- plnorm(x, meanlog, sdlog, lower.tail = FALSE)
    
    #Credible bounds for survival curves
    credint <- res[grepl("pred_comp", row.names(res)), ]
    credint <- credint %>%
      mutate(x = as.numeric(str_extract(row.names(.), "[0-9]+"))) %>%
      select(x, surv_est = `50%`, cilb = `2.5%`, ciub = `97.5%`) %>%
      bind_rows(c(x = 0, surv_est = 1, cilb = 1, ciub = 1))
    
    surv_dens <- cbind.data.frame(x, dens, surv, "label" = label) %>%
      full_join(credint, by = "x")
    
  }else{
    
    #Overall survival and density curves
    sdlog <- param %>% filter(value == "sdlog") %>% pull(est)
    surv_dens <- NULL
    x <- seq(0, 30, 0.1)
    for(sev in c("Minimal", "Moderate", "Advanced")){
      meanlog <- param %>% filter(severity == sev, value == "meanlog") %>% pull(est)
      
      dens <- dlnorm(x, meanlog, sdlog)
      surv <- plnorm(x, meanlog, sdlog, lower.tail = FALSE)
      densTemp <- cbind.data.frame(x, "severity" = sev, meanlog, sdlog, dens, surv)
      
      surv_dens <- bind_rows(surv_dens, densTemp)
    }
    
    #Credible bounds for survival curves
    credint <- res[grepl("pred_[a-z]*", row.names(res)), ]
    credint <- credint %>%
      mutate(x = as.numeric(str_extract(row.names(.), "[0-9]+")),
             severity = ifelse(grepl("min", row.names(.)), "Minimal",
                               ifelse(grepl("mod", row.names(.)), "Moderate",
                                      ifelse(grepl("adv", row.names(.)), "Advanced", NA)))) %>%
      select(x, severity, surv_est = `50%`, cilb = `2.5%`, ciub = `97.5%`) %>%
      bind_rows(cbind.data.frame(x = 0, severity = "Minimal", surv_est = 1, cilb = 1, ciub = 1),
                cbind.data.frame(x = 0, severity = "Moderate", surv_est = 1, cilb = 1, ciub = 1),
                cbind.data.frame(x = 0, severity = "Advanced", surv_est = 1, cilb = 1, ciub = 1))
    
    surv_dens <- surv_dens %>% 
      mutate(label = label) %>%
      full_join(credint, by = c("x", "severity")) %>%
      mutate(severity = factor(severity, levels = c("Minimal", "Moderate", "Advanced", "Unknown"),
                               labels = c("Minimal", "Moderately advanced", "Far advanced", "Unknown")))
  }
  
  return(surv_dens)
}

## Function to format a raw table of study-level results from the Bayesian analysis
format_ind_est <- function(covar, res, data, fixed){
  
  #Individual study median
  med_ind <- as.data.frame(res[grepl("med_ind", row.names(res)),])
  med_ind <- med_ind %>%
    mutate(study_id_num = as.numeric(gsub("med_ind\\[|\\]", "", row.names(med_ind))),
           value = "median")
  
  #Individual study predictions
  pred_ind <- as.data.frame(res[grepl("pred[0-9]*\\[", row.names(res)),])
  pred_ind <- pred_ind %>%
    mutate(rown = row.names(.),
           study_id_num = as.numeric(gsub("pred[0-9]*\\[|\\]", "", rown)),
           value = str_extract(rown, "pred[0-9]*")) %>%
    select(-rown)
  
  #Individual study meanlog
  if(fixed == FALSE){
    mean_ind <- as.data.frame(res[grepl("meanlog", row.names(res)),])
    mean_ind <- mean_ind %>%
      mutate(study_id_num = as.numeric(gsub("meanlog\\[|\\]", "", row.names(mean_ind))),
             value = "meanlog")
  }else{
    mean_ind <- as.data.frame(res[grepl("meanlog_ind", row.names(res)),])
    mean_ind <- mean_ind %>%
      mutate(study_id_num = as.numeric(gsub("meanlog_ind\\[|\\]", "", row.names(mean_ind))),
             value = "meanlog")
  }
  
  #Combining the above
  ind_est <- bind_rows(med_ind, mean_ind, pred_ind)
  names(ind_est) <- c("cilb", "lowerquant", "est", "upperquant", "ciub",
                      "study_id_num", "value")
  ind_est <- ind_est %>%
    full_join(data[[1]], by = "study_id_num") %>%
    left_join(covar, by = "study_id") %>%
    select(-study_id_num)
  
  return(ind_est)
}

## Function to format the study-level survival curves from the Bayesian analysis
format_ind_surv <- function(ind_est, covar, param, label){
  
  #Finding density and survival estimates for each study
  sdlog <- param %>% filter(value == "sdlog") %>% pull(est)
  par_data <- ind_est %>% filter(value == "meanlog")
  
  ind_surv <- NULL
  x <- seq(0, 30, 0.1)
  for(i in 1:nrow(par_data)){
    row <- par_data[i,]
    
    dens <- dlnorm(x, row$est, sdlog)
    surv <- plnorm(x, row$est, sdlog, lower.tail = FALSE)
    densTemp <- cbind.data.frame(x, "study_id" = row$study_id, "meanlog" = row$est,
                                 sdlog, dens, surv)
    
    ind_surv <- bind_rows(ind_surv, densTemp)
  }
  ind_surv <- ind_surv %>%
    left_join(covar, by = "study_id") %>%
    mutate(label = label)
  
  return(ind_surv)
}

## Function to format the study-level predictions (median survival and survival probabilities)
format_pred_comb <- function(ind_est, param, label, fixed){
  
  #Combining study-specific and overall median and predictions
  pred_comb <- ind_est %>% 
    filter(value == "median" | grepl("pred", value)) %>%
    bind_rows(param %>% filter(value == "med" | grepl("pred", value))) %>%
    mutate(shape = ifelse(is.na(study_id), "Overall", "Individual"))
  
  if(fixed == FALSE){
    pred_comb <- pred_comb %>%
      replace_na(list(severity = "",
                      study_id = "Overall",
                      first_author = "Overall"))
  }else{
    pred_comb <- pred_comb %>%
      mutate(study_sev = ifelse(is.na(study_sev), paste("Overall", severity, sep = "_"),
                                study_sev),
             first_author = ifelse(is.na(first_author), "Overall", first_author),
             severity = ifelse(grepl("Overall", study_sev) & severity == "Minimal", "",
                               ifelse(grepl("Overall", study_sev) & severity == "Moderate", " ",
                                      ifelse(grepl("Overall", study_sev) & severity == "Advanced", "  ",
                                             severity))),
             severity = factor(severity, levels = c("Minimal", "", "Moderate", " ", "Advanced", "  "),
                               labels = c("Minimal", "", "Moderately\nadvanced", " ", "Far\nadvanced", "  ")))
    
  }
  
  pred_comb <- pred_comb %>% 
    mutate(pred_label = factor(value, levels = c("pred1", "pred5", "pred10", "median"),
                               labels = c("1-Year", "5-Year", "10-Year", "Median")),
           label = label)
  
  return(pred_comb)
}


## Function to format results of Bayesian analysis
formatBayesian <- function(mortalityData, res, data, label, fixed = FALSE){
  
  #Parameter table
  param <- format_param(res, label, fixed)
  
  #Overall survival and density table
  surv_dens <- format_surv_dens(param, res, label, fixed)
  
  #Creating table of covariate data
  covar <- mortalityData %>%
    #full_join(studyid, by = "study_id") %>%
    group_by(study_id) %>%
    summarize(first_author = first(first_author),
              sanatorium = first(sanatorium),
              location = first(location),
              time_period = first(time_period),
              start_type = first(start_type),
              .groups = "drop")
  
  #Study-level raw results
  ind_est <- format_ind_est(covar, res, data, fixed)

  #Study-level survival curves
  ind_surv <- format_ind_surv(ind_est, covar, param, label)

  #Study-level predictions (median survival and survival probabilities)
  pred_comb <- format_pred_comb(ind_est, param, label, fixed)
  
  return(list("surv_dens" = surv_dens, "param" = param,
              "pred_comb" = pred_comb, "ind_surv" = ind_surv))
}


## Function to format cure analysis tables
make_cure_tab <- function(eval_tab, aggregate_tab){
  
  #Odds ratios for minimal and moderate vs. advanced
  or <- as.data.frame(summary(eval_tab[, c("ORmin", "ORmod")])$quantiles) %>%
    mutate(rownames = row.names(.),
           severity = ifelse(grepl("min", rownames), "Minimal", "Moderate"),
           OR_CI = paste0(round(`50%`, 2), " (", round(`2.5%`, 2), ", ", round(`97.5%`, 2), ")")) %>%
    select(severity, OR_CI)
  
  #Table of counts per study
  cureTab <- aggregate_tab %>%
    mutate(study_id = as.character(study_id)) %>%
    left_join(studyid, by = "study_id") %>%
    mutate(pMin = 100 * round(cMin / nMin, 2),
           pMod = 100 * round(cMod / nMod, 2),
           pAdv = 100 * round(cAdv / nAdv, 2),
           Min_cure = ifelse(nMin == 0, "-", paste0(cMin, " (", pMin, "%)")),
           Mod_cure = paste0(cMod, " (", pMod, "%)"),
           Adv_cure = paste0(cAdv, " (", pAdv, "%)")) %>%
    select(first_author, Min_total = nMin, Min_cure, Mod_total = nMod, Mod_cure,
           Adv_total = nAdv, Adv_cure) %>%
    arrange(first_author)
  
  return(list(or, cureTab))
}


dblTochr <- function(dataset){
  dataset$study_id <- as.character(dataset$study_id)
  dataset$paper_id <- as.character(dataset$paper_id)
  
  return(dataset)
}

calcCureRate <- function(cureData){
  cureData$cureRate <- cureData$c2/cureData$n
  
  return(cureData)
  
}


