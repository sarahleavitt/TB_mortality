#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program contains the functions to perform the mortality Bayesian 
# meta-analysis with a fixed effect for strata in separate models.
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


#### MCMC Model Function -------------------------------------------------------

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
  #Overall median survival for one unit change in strata value
  med_timep <- exp(alpha + btimep)
  med_loc <- exp(alpha + bloc)
  med_san <- exp(alpha + bsan)
  
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
                    "med_timep", "med_loc", "med_san",
                    "meanlog_ind", "med_ind","pred1", "pred5", "pred10")


#### Function to Run the Models ------------------------------------------------

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


#### Running the Model ---------------------------------------------------------

#Creating interval variables
#interval = 1 implies x1 < t <= x2
#dinterval() isn't working as it should so make all obs have interval = 1
#and set x2 to be 10000 (close enough to infinity) for right censored
mortality <- mortality %>% 
  mutate(time = ifelse(death_tb == 0, NA, interval_l),
         interval = 1,
         x1 = interval_l,
         x2 = ifelse(death_tb == 0, 10000, interval_r))

# n.iter <- 31000
# n.burnin <- 1000
# n.thin <- 30

n.iter <- 100
n.burnin <- 10
n.thin <- 1


#Subsetting and formatting data
mortality_strata <- mortality %>%
  #Removing severity stratified mortality data for study 79_1023 (only 4-year follow-up)
  filter(!cohort_id %in% c("79_1023_5", "79_1023_6", "79_1023_7")) %>%
  mutate(study_id_num = as.numeric(factor(study_id)),
         pre1930 = as.numeric(time_period == "pre-1930"),
         northamerica = as.numeric(location == "North America"),
         sanatorium = as.numeric(sanatorium == "Yes"))

#Getting information for each study
study_data <- mortality_strata %>%
  group_by(study_id_num) %>%
  summarize(pre1930 = first(pre1930),
            northamerica = first(northamerica),
            sanatorium = first(sanatorium))

output_all <- run_all_strata(mortality_strata, study_data, n.iter = n.iter,
                             n.burnin = n.burnin, n.thin = n.thin)


#### Formatting and Saving Results ---------------------------------------------

data <- getData(mortality_strata)
res_all <- output_all$res
eval_all <- output_all$eval

save(data, res_all, eval_all, file = "R/bayesian_allstrata.RData")

png("Figures/xyplot_all.png")
xyplot(eval_all)
dev.off()
png("Figures/autocorr_all.png")
autocorr.plot(eval_all)
dev.off()



