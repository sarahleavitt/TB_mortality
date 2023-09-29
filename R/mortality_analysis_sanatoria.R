# Rerun analysis on sanatoria stratified by start of follow-up: entry or exit 

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

#Subsetting and formatting data
mortality_strata <- mortality %>%
  #Removing severity stratified mortality data for study 79_1023 (only 4-year follow-up)
  filter(!cohort_id %in% c("79_1023_5", "79_1023_6", "79_1023_7")) %>%
  mutate(study_id_num = as.numeric(factor(study_id)))

#### Running the Models --------------------------------------------------------
## Sanatorium at entry and sanitorium starting at exit
### make sure to place Ferguson as exit study since only those who survived were tracked
san_entry <- mortality_strata %>%
  filter(sanatorium == "Yes" & start_type == "Entry" & study_id != "94") %>%
  mutate(study_id_num = as.numeric(factor(study_id)))

san_exit <- mortality_strata %>%
  filter((sanatorium == "Yes" & start_type == "Exit") | study_id == "94") %>%
  mutate(study_id_num = as.numeric(factor(study_id)))

output_san_entry <- run_comp(san_entry, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
output_san_exit <- run_comp(san_exit, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)

#### Formatting and Saving Results ---------------------------------------------

data_san_entry <- getData(san_entry)
res_san_entry <- output_san_entry$res
eval_san_entry <- output_san_entry$eval
data_san_exit <- getData(san_exit)
res_san_exit <- output_san_exit$res
eval_san_exit <- output_san_exit$eval

save(data_san_entry, res_san_entry, eval_san_entry,
     data_san_exit, res_san_exit, eval_san_exit, file = "R/bayesian_raw_entry.RData")

##############################################################################
# This program formats the mortality analysis results and saves the
# diagnostic plots
##############################################################################

#Reading in the study_id correspondence table
studyid <- read.csv("data/study_id.csv")

## Formatting results
form_san_entry <- formatBayesian(mortality, res_san_entry, data_san_entry, "Sanatorium entry")
form_san_exit <- formatBayesian(mortality, res_san_exit, data_san_exit, "Sanatorium exit")

## Saving results
save(form_san_entry,form_san_exit, file = "R/bayesian_clean_sanitoria.RData")

#### Diagnostic Plots ------------------------------------------------------------------------------

#Entry analysis
png("Figures/xyplot_san_entry.png")
xyplot(eval_san_entry)
dev.off()
png("Figures/autocorr_san_entry.png")
autocorr.plot(eval_san_entry)
dev.off()

#Exit analysis
png("Figures/xyplot_san_exit.png")
xyplot(eval_san_exit)
dev.off()
png("Figures/autocorr_san_exit.png")
autocorr.plot(eval_san_exit)
dev.off()

####################################################
# Create figures of the mortality analysis results
#####################################################

#Reading in analysis results
#load('R/bayesian_clean_sanitoria.RData')

#### Survival Curves: Stratified -----------------------------------------------

# Combined plot for all stratifications
surv_all <- bind_rows(form_san_entry$surv_dens, form_san_exit$surv_dens) %>%
  mutate(label = factor(label, levels = c("Sanatorium entry",
                                          "Sanatorium exit")))

ind_all <- bind_rows(form_san_entry$ind_surv, form_san_exit$ind_surv) %>%
  mutate(label = factor(label, levels = c("Sanatorium entry",
                                          "Sanatorium exit")))
ggplot(surv_all) +
  geom_line(aes(x = x, y = surv),
            color = "black", size = 1, linetype = "longdash") +
  geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
              stat = "identity", linetype = 0, alpha = 0.4) +
  facet_wrap(~label, nrow = 2) +
  geom_line(data = ind_all,
            aes(x = x, y = surv, group = study_id),
            size = 0.6, alpha = 0.2) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw()

ggsave("Figures/survival_curve_sanitoria.png", width = 6.5, height = 5)

##############################################################################
# Create tables of the mortality analysis results
##############################################################################
#load('R/bayesian_raw.RData')
#load('R/bayesian_clean.RData')

#### Table of Main Results -----------------------------------------------------

#Combining raw tables from results lists
raw_tab <- bind_rows(form_san_entry$param,
                     form_san_exit$param)

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

data_san_entry2 <- as.data.frame(t(data_san_entry[[2]]))
data_san_entry2$label <- "Sanatorium entry"
data_san_exit2 <- as.data.frame(t(data_san_exit[[2]]))
data_san_exit2$label <- "Sanatorium exit"

counts <- bind_rows(data_san_entry2, data_san_exit2) %>%
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

write.csv(final_tab, "data/mortality_results_table_sanitoria.csv")

