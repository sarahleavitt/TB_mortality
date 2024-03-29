#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Mortality Analysis

##############################################################################
# This program contains the performs the mortality Bayesian meta-analysis 
# for the combined model and stratified models.
# This program takes about 1 hour to run.
##############################################################################

options(scipen=999)
options(digits = 10)
set.seed(150183)

rm(list = ls())
source("R/utils.R")
reload_source()

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
  #Removing people who are censored at the study start
  filter(no_survival == FALSE) %>%
  mutate(time = ifelse(death_tb == 0, NA, interval_l),
         interval = 1,
         x1 = interval_l,
         x2 = ifelse(death_tb == 0, 10000, interval_r))

n.iter <- 31000
n.burnin <- 1000
n.thin <- 30

# n.iter <- 100
# n.burnin <- 10
# n.thin <- 1

#Subsetting and formatting data
mortality_strata <- mortality %>%
  #Removing severity stratified mortality data for study 79_1023 (only 4-year follow-up)
  filter(!cohort_id %in% c("79_1023_5", "79_1023_6", "79_1023_7")) %>%
  mutate(study_id_num = as.numeric(factor(study_id)))


#### Running the Models --------------------------------------------------------

## Combined Model
output_comb <- run_comp(mortality_strata, n.iter = n.iter, n.burnin = n.burnin,
                        n.thin = n.thin)


## Time period stratification
pre <- mortality_strata %>%
  filter(time_period == "pre-1930") %>%
  mutate(study_id_num = as.numeric(factor(study_id)))

post <- mortality_strata %>%
  filter(time_period == "post-1930") %>%
  mutate(study_id_num = as.numeric(factor(study_id)))

output_pre <- run_comp(pre, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
output_post <- run_comp(post, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)


## Location stratification
namerica <- mortality_strata %>%
  filter(location == "North America") %>%
  mutate(study_id_num = as.numeric(factor(study_id)))

europe <- mortality_strata %>%
  filter(location == "Europe") %>%
  mutate(study_id_num = as.numeric(factor(study_id)))

output_namerica <- run_comp(namerica, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
output_europe <- run_comp(europe, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)


## Sanatorium stratification
yessan <- mortality_strata %>%
  filter(sanatorium == "Yes") %>%
  mutate(study_id_num = as.numeric(factor(study_id)))

nosan <- mortality_strata %>%
  filter(sanatorium == "No") %>%
  mutate(study_id_num = as.numeric(factor(study_id)))

output_yessan <- run_comp(yessan, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
output_nosan <- run_comp(nosan, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)


#### Formatting and Saving Results ---------------------------------------------

data_comb <- getData(mortality_strata)
res_comb <- output_comb$res
eval_comb <- output_comb$eval

data_pre <- getData(pre)
res_pre <- output_pre$res
eval_pre <- output_pre$eval
data_post <- getData(post)
res_post <- output_post$res
eval_post <- output_post$eval

data_namerica <- getData(namerica)
res_namerica <- output_namerica$res
eval_namerica <- output_namerica$eval
data_europe <- getData(europe)
res_europe <- output_europe$res
eval_europe <- output_europe$eval

data_yessan <- getData(yessan)
res_yessan <- output_yessan$res
eval_yessan <- output_yessan$eval
data_nosan <- getData(nosan)
res_nosan <- output_nosan$res
eval_nosan <- output_nosan$eval

save(data_comb, res_comb, eval_comb,
     data_pre, res_pre, eval_pre,
     data_post, res_post, eval_post,
     data_namerica, res_namerica, eval_namerica,
     data_europe, res_europe, eval_europe,
     data_yessan, res_yessan, eval_yessan,
     data_nosan, res_nosan, eval_nosan, file = "R/bayesian_raw.RData")
