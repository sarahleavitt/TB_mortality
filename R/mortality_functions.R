#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program contains the functions to perform the mortality Bayesian 
# meta-analysis.
##############################################################################



#### MCMC Model Functions --------------------------------------------------------------------------

#### Complete model (no fixed effect for severity) ####

m_comp <- function(){
  
  for(j in 1:n_frail){
    
    #Frailty for each study
    ui[j] ~ dnorm(mu, tau)
    #Meanlog of lognormal distribution for each study
    meanlog[j] <- ui[j]
    #Median of lognormal distribution for each study
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
  #sdlog of all survival densities (study-specific and overall)
  sdlog <- sqrt(1/taulog)
  #Overall median
  med_comp <- exp(mu)
}

#Parameters to track
par_comp <- c("mu", "theta", "sdlog", "med_comp", "meanlog", "med_ind",
              "pred_comp", "pred1", "pred5", "pred10")


#### Stratified model (fixed effect for severity) ####

m_sev <- function(){
  
  #Frailty for each study
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
  #sdlog of all survival densities (study-specific and overall)
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
  
  #study-specific meanlog of the lognormal density and median survival
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

run_comp <- function(df, n.iter = 31000, n.burnin = 1000, n.thin = 30){
  
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
  
  #Getting information for each study
  study_data <- df %>%
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
             n_study_sev = nrow(study_data),
             frail2 = study_data$study_id_num,
             study_sev_mod = study_data$sev_mod,
             study_sev_adv = study_data$sev_adv
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



#### Formatting results ----------------------------------------------------------------------------

## Function to get information about the dataset for each run
getData <- function(data){
  
  #Finding concordance to frailty IDs
  tab <- data %>%
    select(study_id, study_id_num, study_sev, study_sev_num) %>%
    filter(!duplicated(study_sev))
  
  #Finding n studies, severity groups, studys, individuals
  counts <- c("nStudies" = length(unique(data$study_id_num)),
              "nSeverity" = length(unique(data$study_sev)),
              "nCohorts" = length(unique(data$cohort_id)),
              "nIndividuals" = nrow(data))
  
  return(list(tab, counts))
}

