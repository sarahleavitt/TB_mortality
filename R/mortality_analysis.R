#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program performs the mortality Bayesian meta-analysis and saves the
# results as an R object. This takes around 3 hours to run.
##############################################################################

options(scipen=999)
options(digits = 10)
set.seed(150183)

rm(list = ls())
source("R/utils.R")
reload_source()

#Reading in individual mortality data
mortality <- read.csv("data/mortality_data.csv")

#Creating interval variables
#interval = 1 implies x1 < t <= x2
#dinterval() isn't working as it should so make all obs have interval = 1
#and set x2 to be 10000 (close enough to infinity) for right censored
mortality <- mortality %>% 
  mutate(time = ifelse(death_tb == 0, NA, interval_l),
         interval = 1,
         x1 = interval_l,
         x2 = ifelse(death_tb == 0, 10000, interval_r))



#### MCMC Model Functions --------------------------------------------------------------------------

#### Complete model (no fixed effect for severity) ####

m_comp <- function(){
  
  for(j in 1:n_frail){
    
    #Frailty for each cohort
    ui[j] ~ dnorm(mu, tau)
    #Meanlog of lognormal distribution for each cohort
    meanlog[j] <- ui[j]
    #Median of lognormal distribution for each cohort
    med_ind[j] <- exp(meanlog[j])
    
    #Prediction of 1/5/10 year survival
    pred1[j] <- 1 - plnorm(1, meanlog[j], taulog)
    pred5[j] <- 1 - plnorm(5, meanlog[j], taulog)
    pred10[j] <- 1 - plnorm(10, meanlog[j], taulog)
  }
  
  for(i in 1:N){
    #Setting up the survival model
    interval[i] ~ dinterval(time[i], lim[i, ])
    time[i] ~ dlnorm(meanlog[frail[i]], taulog)
  }
  
  #Priors
  mu ~ dnorm(0, 0.0001)
  taulog ~ dgamma(1, 1)
  tau ~ dgamma(1, 1)
  
  #Prediction for confidence bands
  for(t in 1:30){
    pred_comp[t] <- 1 - plnorm(t, mu, taulog)
  }
  
  ## Derived parameters ##
  
  #Variance of the frailty distribution
  theta <- 1/tau
  #sdlog of all survival densities (cohort-specific and overall)
  sdlog <- sqrt(1/taulog)
  #Overall median
  med_comp <- exp(mu)
}

#Parameters to track
par_comp <- c("mu", "theta", "sdlog", "med_comp", "meanlog", "med_ind",
             "pred_comp", "pred1", "pred5", "pred10")



#### Stratified model (fixed effect for severity) ####

m_sev <- function(){
  
  #Frailty for each cohort
  for(j in 1:n_frail){
    ui[j] ~ dnorm(alpha, tau)
  }
  
  for(i in 1:N){
    #Meanlog of lognormal distribution for each individual
    meanlog[i] <- ui[frail[i]] + bmod*sev_mod[i] + badv*sev_adv[i]
    #Setting up the survival model
    interval[i] ~ dinterval(time[i], lim[i, ])
    time[i] ~ dlnorm(meanlog[i], taulog)
  }
  
  #Priors
  alpha ~ dnorm(0, 0.0001)
  bmod ~ dnorm(0, 0.0001)
  badv ~ dnorm(0, 0.0001)
  taulog ~ dgamma(1, 1)
  tau ~ dgamma(1, 1)
  
  #Variance of the frailty distribution
  theta <- 1/tau
  #sdlog of all survival densities (cohort-specific and overall)
  sdlog <- sqrt(1/taulog)
  #Overall median survival for each severity
  meanlog_min <- alpha
  meanlog_mod <- alpha + bmod
  meanlog_adv <- alpha + badv
  med_min <- exp(alpha)
  med_mod <- exp(alpha + bmod)
  med_adv <- exp(alpha + badv)
  
  #Prediction for confidence bands
  for(t in 1:30){
    pred_min[t] <- 1 - plnorm(t, meanlog_min, taulog)
    pred_mod[t] <- 1 - plnorm(t, meanlog_mod, taulog)
    pred_adv[t] <- 1 - plnorm(t, meanlog_adv, taulog)
  }
  
  #Cohort-specific meanlog of the lognormal density and median survival
  for(k in 1:n_study_sev){
    meanlog_ind[k] <- ui[frail2[k]] + bmod*study_sev_mod[k] + badv*study_sev_adv[k]
    med_ind[k] <- exp(meanlog_ind[k])
    
    #Prediction of 1/5/10 year survival
    pred1[k] <- 1 - plnorm(1, meanlog_ind[k], taulog)
    pred5[k] <- 1 - plnorm(5, meanlog_ind[k], taulog)
    pred10[k] <- 1 - plnorm(10, meanlog_ind[k], taulog)
  }
}

#Parameters to track
par_sev <- c("theta", "sdlog", "alpha", "bmod", "badv",
             "meanlog_ind", "meanlog_min", "meanlog_mod", "meanlog_adv",
             "med_ind", "med_min", "med_mod", "med_adv",
             "pred1", "pred5", "pred10", "pred_min", "pred_mod", "pred_adv")




#### Functions to Run the Models -------------------------------------------------------------------


#### Complete model (no fixed effect for severity) ####

run_comp <- function(df, n.iter = 61000, n.burnin = 1000, n.thin = 30){
  
  #Create MCMC dataset
  dt <- list(N = nrow(df),
             interval = df$interval,
             lim = cbind(df$x1, df$x2),
             time = rep(NA, nrow(df)),
             n_frail = length(unique(df$study_sev_num)),
             frail = df$study_sev_num
  )
  
  #Fitting the model
  fit <- jags(data = dt, model.file = m_comp,
              parameters.to.save = par_comp,
              n.iter = n.iter, n.burnin = n.burnin,
              n.chains = 1, n.thin = n.thin)
  
  #Extracting results
  mcmc <- as.mcmc(fit)
  eval <- mcmc[, c("mu", "theta", "sdlog")]
  res <- as.data.frame(summary(mcmc)$quantiles)
  
  return(list("res" = res, "eval" = eval))
}


#### Stratified model (fixed effect for severity) ####

run_sev <- function(df, n.iter = 61000, n.burnin = 1000, n.thin = 30){
  
  #Getting information for each cohort
  cohort_data <- df %>%
    group_by(study_sev_num) %>%
    summarize(study_id_num = first(study_id_num),
              sev_mod = first(sev_mod),
              sev_adv = first(sev_adv),
              sev_unk = first(sev_unk))
  
  #Creating MCMC dataset
  dt <- list(N = nrow(df),
             interval = df$interval,
             lim = cbind(df$x1, df$x2),
             time = rep(NA, nrow(df)),
             n_frail = length(unique(df$study_id_num)),
             frail = df$study_id_num,
             sev_mod = df$sev_mod,
             sev_adv = df$sev_adv,
             n_study_sev = nrow(cohort_data),
             frail2 = cohort_data$study_id_num,
             study_sev_mod = cohort_data$sev_mod,
             study_sev_adv = cohort_data$sev_adv
  )
  
  #Fitting model
  fit <- jags(data = dt, model.file = m_sev,
              parameters.to.save = par_sev,
              n.iter = n.iter, n.burnin = n.burnin,
              n.chains = 1, n.thin = n.thin)
  
  #Extracting results
  mcmc <- as.mcmc(fit)
  eval <- mcmc[, c("alpha", "bmod", "badv", "theta", "sdlog")]
  res <- as.data.frame(summary(mcmc)$quantiles)
  
  return(list("res" = res, "eval" = eval))
}



#### All Studies -----------------------------------------------------------------------------------


#### All Studies: Complete Model ####

#Subsetting and formatting data
mortality_comp_all <- mortality %>%
  #Removing severity stratified mortality data for study 79_1023 (only 4-year follow-up)
  filter(!cohort_id %in% c("79_1023_5", "79_1023_6", "79_1023_7")) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)))

#Running model
output_comp_all <- run_comp(mortality_comp_all)

#Saving diagnostic plots
png("Figures/xyplot_comp_all.png")
xyplot(output_comp_all$eval)
dev.off()
png("Figures/autocorr_comp_all.png")
autocorr.plot(output_comp_all$eval)
dev.off()



#### All Studies: Stratified Model ####

#Subsetting and formatting data
mortality_sev_all <- mortality %>%
  filter(severity != "Unknown") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)),
         sev_mod = as.numeric(severity == "Moderate"),
         sev_adv = as.numeric(severity == "Advanced"),
         sev_unk = as.numeric(severity == "Unknown"))

#Running model
output_sev_all <- run_sev(mortality_sev_all)

#Saving diagnostic plots
png("Figures/xyplot_sev_all.png")
xyplot(output_sev_all$eval)
dev.off()
png("Figures/autocorr_sev_all.png")
autocorr.plot(output_sev_all$eval)
dev.off()



#### US Studies-------------------------------------------------------------------------------------


#### US Studies: Complete Model ####

#Subsetting and formatting data
mortality_comp_us <- mortality %>%
  filter(study_id %in% c("1029", "93", "45", "63", "67", "90_1016")) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)))

#Running model
output_comp_us <- run_comp(mortality_comp_us)

#Saving diagnostic plots
png("Figures/xyplot_comp_us.png")
xyplot(output_comp_us$eval)
dev.off()
png("Figures/autocorr_comp_us.png")
autocorr.plot(output_comp_us$eval)
dev.off()



#### US Studies: Stratified Model ####

#Subsetting and formatting data
mortality_sev_us <- mortality %>%
  filter(study_id %in% c("1029", "93", "45", "63", "67", "90_1016")) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)),
         sev_mod = as.numeric(severity == "Moderate"),
         sev_adv = as.numeric(severity == "Advanced"),
         sev_unk = as.numeric(severity == "Unknown"))

#Running model
output_sev_us <- run_sev(mortality_sev_us)

#Saving diagnostic plots
png("Figures/xyplot_sev_us.png")
xyplot(output_sev_us$eval)
dev.off()
png("Figures/autocorr_sev_us.png")
autocorr.plot(output_sev_us$eval)
dev.off()



#### US post-1930 Studies---------------------------------------------------------------------------


#### US post-1930 Studies: Complete Model ####

#Removing severity stratified but only four-year mortality data for study 79_1023
mortality_comp_post <- mortality %>%
  filter(study_id %in% c("1029", "93", "45")) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)))

#Running model
output_comp_post <- run_comp(mortality_comp_post)

#Saving diagnostic plots
png("Figures/xyplot_comp_post.png")
xyplot(output_comp_post$eval)
dev.off()
png("Figures/autocorr_comp_post.png")
autocorr.plot(output_comp_post$eval)
dev.off()



#### US post-1930 Studies: Stratified Model ####

#Subsetting and formatting data
mortality_sev_post <- mortality %>%
  filter(study_id %in% c("1029", "93", "45")) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)),
         sev_mod = as.numeric(severity == "Moderate"),
         sev_adv = as.numeric(severity == "Advanced"),
         sev_unk = as.numeric(severity == "Unknown"))

#Running model
output_sev_post <- run_sev(mortality_sev_post)

#Saving diagnostic plots
png("Figures/xyplot_sev_post.png")
xyplot(output_sev_post$eval)
dev.off()
png("Figures/autocorr_sev_post.png")
autocorr.plot(output_sev_post$eval)
dev.off()



#### Sanatorium Sensitivity Analysis----------------------------------------------------------------


#### Sanatorium/hospital ####

#Subsetting data
san <- mortality_comp_all %>%
  filter(sanatorium == "Yes") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)))

#Running model
output_san <- run_comp(san)

#Making diagnostic plots
png("Figures/xyplot_san.png")
xyplot(output_san$eval)
dev.off()
png("Figures/autocorr_san.png")
autocorr.plot(output_san$eval)
dev.off()



#### TB mortality: non-sanatorium ####

#Subsetting data
nosan <- mortality_comp_all %>%
  filter(sanatorium == "No") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)))

#Running model
output_nosan <- run_comp(nosan)

#Making diagnostic plots
png("Figures/xyplot_nosan.png")
xyplot(output_nosan$eval)
dev.off()
png("Figures/autocorr_nosan.png")
autocorr.plot(output_nosan$eval)
dev.off()




#### Compiling Results------------------------------------------------------------------------------


## Function to get information about the dataset for each run
getData <- function(data){

  #Finding concordance to frailty IDs
  tab <- data %>%
    select(study_id, study_id_num, study_sev, study_sev_num) %>%
    filter(!duplicated(study_sev))
  
  #Finding n studies, severity groups, cohorts, individuals
  counts <- c("nStudies" = length(unique(data$study_id_num)),
                 "nSeverity" = length(unique(data$study_sev)),
                 "nCohorts" = length(unique(data$cohort_id)),
                 "nIndividuals" = nrow(data))
  
  return(list(tab, counts))
}

data_comp_all <- getData(mortality_comp_all)
data_sev_all <- getData(mortality_sev_all)
data_comp_us  <- getData(mortality_comp_us)
data_sev_us <- getData(mortality_sev_us)
data_comp_post  <- getData(mortality_comp_post)
data_sev_post <- getData(mortality_sev_post)
data_san <- getData(san)
data_nosan <- getData(nosan)

res_comp_all <- output_comp_all$res
res_sev_all <- output_sev_all$res
res_comp_us <- output_comp_us$res
res_sev_us <- output_sev_us$res
res_comp_post <- output_comp_post$res
res_sev_post <- output_sev_post$res
res_san <- output_san$res
res_nosan <- output_nosan$res

eval_comp_all <- output_comp_all$eval
eval_sev_all <- output_sev_all$eval
eval_comp_us <- output_comp_us$eval
eval_sev_us <- output_sev_us$eval
eval_comp_post <- output_comp_post$eval
eval_sev_post <- output_sev_post$eval
eval_san <- output_san$eval
eval_nosan <- output_nosan$eval


save(res_comp_all, res_sev_all, eval_comp_all, eval_sev_all, data_comp_all, data_sev_all,
     res_comp_us, res_sev_us, eval_comp_us, eval_sev_us, data_comp_us, data_sev_us,
     res_comp_post, res_sev_post, eval_comp_post, eval_sev_post, data_comp_post, data_sev_post,
     res_san, res_nosan, eval_san, eval_nosan, data_san, data_nosan,
     file = "R/bayesian_mortality.RData")

