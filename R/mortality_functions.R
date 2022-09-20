#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program contains the functions to perform the mortality Bayesian 
# meta-analysis.
##############################################################################



#### MCMC Model Functions --------------------------------------------------------------------------

#### Complete model (no fixed effects for strata) ####

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


#### Fixed effect for strata separately ####

m_one_strata <- function(){
  
  #Frailty for each study
  for(j in 1:n_frail){
    ui[j] ~ dnorm(alpha, tau)
  }
  
  for(i in 1:N){
    #Meanlog of lognormal distribution for each individual
    meanlog[i] <- ui[frail[i]] + bstrata*strata[i]
    #Setting up the survival model
    interval[i] ~ dinterval(time[i], lim[i, ])
    time[i] ~ dlnorm(meanlog[i], taulog)
  }
  
  #Priors
  alpha ~ dnorm(0, 0.0001)
  bstrata ~ dnorm(0, 0.0001)
  taulog ~ dgamma(1, 1)
  tau ~ dgamma(1, 1)
  
  #Variance of the frailty distribution
  theta <- 1/tau
  #sdlog of all survival densities (study-specific and overall)
  sdlog <- sqrt(1/taulog)
  #Overall median survival for each severity
  meanlog_strata0 <- alpha
  meanlog_strata1 <- alpha + bstrata
  med_strata0 <- exp(alpha)
  med_strata1 <- exp(alpha + bstrata)
  
  #Prediction for confidence bands
  for(t in 1:30){
    pred_strata0[t] <- 1 - plnorm(t, meanlog_strata0, taulog)
    pred_strata1[t] <- 1 - plnorm(t, meanlog_strata1, taulog)
  }
  
  #study-specific meanlog of the lognormal density and median survival
  for(k in 1:n_study){
    meanlog_ind[k] <- ui[frail2[k]] + bstrata*study_strata[k]
    med_ind[k] <- exp(meanlog_ind[k])
    
    #Prediction of 1/5/10 year survival
    pred1[k] <- 1 - plnorm(1, meanlog_ind[k], taulog)
    pred5[k] <- 1 - plnorm(5, meanlog_ind[k], taulog)
    pred10[k] <- 1 - plnorm(10, meanlog_ind[k], taulog)
  }
}

#Parameters to track
par_one_strata <- c("theta", "sdlog", "alpha", "bstrata",
                    "meanlog_ind", "meanlog_strata0", "meanlog_strata1",
                    "med_ind", "med_strata0", "med_strata1",
                    "pred1", "pred5", "pred10", "pred_strata0", "pred_strata1")



#### Fixed effect for strata combined ####

m_all_strata <- function(){
  
  #Frailty for each study
  for(j in 1:n_frail){
    ui[j] ~ dnorm(alpha, tau)
  }
  
  for(i in 1:N){
    #Meanlog of lognormal distribution for each individual
    meanlog[i] <- ui[frail[i]] + btimep*timep[i] + bloc*loc[i] + bsan*san[i]
    #Setting up the survival model
    interval[i] ~ dinterval(time[i], lim[i, ])
    time[i] ~ dlnorm(meanlog[i], taulog)
  }
  
  #Priors
  alpha ~ dnorm(0, 0.0001)
  btimep ~ dnorm(0, 0.0001)
  bloc ~ dnorm(0, 0.0001)
  bsan ~ dnorm(0, 0.0001)
  taulog ~ dgamma(1, 1)
  tau ~ dgamma(1, 1)
  
  #Variance of the frailty distribution
  theta <- 1/tau
  #sdlog of all survival densities (study-specific and overall)
  sdlog <- sqrt(1/taulog)
  
  #study-specific meanlog of the lognormal density and median survival
  for(k in 1:n_study){
    meanlog_ind[k] <- ui[frail2[k]] + btimep*study_timep[k] + bloc*study_loc[k] + bsan*study_san[k]
    med_ind[k] <- exp(meanlog_ind[k])
    
    #Prediction of 1/5/10 year survival
    pred1[k] <- 1 - plnorm(1, meanlog_ind[k], taulog)
    pred5[k] <- 1 - plnorm(5, meanlog_ind[k], taulog)
    pred10[k] <- 1 - plnorm(10, meanlog_ind[k], taulog)
  }
}

#Parameters to track
par_all_strata <- c("theta", "sdlog", "alpha", "btimep", "bloc", "bsan",
                    "meanlog_ind", "med_ind","pred1", "pred5", "pred10")




#### Functions to Run the Models -------------------------------------------------------------------


#### Complete model (no fixed effects for strata) ####

run_comp <- function(df, n.iter = 31000, n.burnin = 1000, n.thin = 30){
  
  #Create MCMC dataset
  dt <- list(N = nrow(df),
             interval = df$interval,
             lim = cbind(df$x1, df$x2),
             time = rep(NA, nrow(df)),
             n_frail = length(unique(df$study_id_num)),
             frail = df$study_id_num
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


#### Fixed effect for strata separately ####

run_one_strata <- function(df, study_data, strata_var, n.iter = 61000,
                           n.burnin = 1000, n.thin = 30){
  
  study_data <- as.data.frame(study_data)
  
  #Creating MCMC dataset
  dt <- list(N = nrow(df),
             interval = df$interval,
             lim = cbind(df$x1, df$x2),
             time = rep(NA, nrow(df)),
             n_frail = length(unique(df$study_id_num)),
             frail = df$study_id_num,
             strata = df[, strata_var],
             n_study = nrow(study_data),
             frail2 = study_data$study_id_num,
             study_strata = study_data[, strata_var]
  )
  
  #Fitting model
  fit <- jags(data = dt, model.file = m_one_strata,
              parameters.to.save = par_one_strata,
              n.iter = n.iter, n.burnin = n.burnin,
              n.chains = 1, n.thin = n.thin)
  
  #Extracting results
  mcmc <- as.mcmc(fit)
  eval <- mcmc[, c("alpha", "bstrata", "theta", "sdlog")]
  res <- as.data.frame(summary(mcmc)$quantiles)
  
  return(list("res" = res, "eval" = eval))
}


#### Fixed effect for strata combined ####

run_all_strata <- function(df, study_data, n.iter = 61000, n.burnin = 1000,
                           n.thin = 30){
  
  #Creating MCMC dataset
  dt <- list(N = nrow(df),
             interval = df$interval,
             lim = cbind(df$x1, df$x2),
             time = rep(NA, nrow(df)),
             n_frail = length(unique(df$study_id_num)),
             frail = df$study_id_num,
             timep = df$pre1930,
             loc = df$northamerica,
             san = df$sanatorium,
             n_study = nrow(study_data),
             frail2 = study_data$study_id_num,
             study_timep = study_data$pre1930,
             study_loc = study_data$northamerica,
             study_san = study_data$sanatorium
  )
  
  #Fitting model
  fit <- jags(data = dt, model.file = m_all_strata,
              parameters.to.save = par_all_strata,
              n.iter = n.iter, n.burnin = n.burnin,
              n.chains = 1, n.thin = n.thin)
  
  #Extracting results
  mcmc <- as.mcmc(fit)
  eval <- mcmc[, c("alpha", "btimep", "bloc", "bsan", "theta", "sdlog")]
  res <- as.data.frame(summary(mcmc)$quantiles)
  
  return(list("res" = res, "eval" = eval))
}



#### Formatting results ----------------------------------------------------------------------------

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
