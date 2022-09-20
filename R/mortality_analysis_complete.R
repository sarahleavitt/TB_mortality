#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program contains the functions to perform the mortality Bayesian 
# meta-analysis for the complete model (no stratification).
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



#### Function to Run the Models ------------------------------------------------

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
  mutate(study_id_num = as.numeric(factor(study_id)))

#Running the model
output_comp <- run_comp(mortality_strata, n.iter = n.iter, n.burnin = n.burnin,
                        n.thin = n.thin)


#### Formatting and Saving Results ---------------------------------------------

data <- getData(mortality_strata)
res_comp <- output_comp$res
eval_comp <- output_comp$eval

save(data, res_comp, eval_comp, file = "R/bayesian_comp.RData")

form_comp <- formatBayesian(mortality, res_comp, data, "Combined")

png("Figures/xyplot_comp.png")
xyplot(eval_comp)
dev.off()
png("Figures/autocorr_comp.png")
autocorr.plot(eval_comp)
dev.off()

