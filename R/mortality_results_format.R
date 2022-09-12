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

#Reading in individual mortality data
mortality <- read.csv("data/mortality_data.csv")

#Reading in analysis results
load('R/bayesian_all.RData')
load('R/bayesian_us.RData')
load('R/bayesian_nonus.RData')
load('R/bayesian_post.RData')
load('R/bayesian_pre.RData')
load('R/bayesian_san.RData')

## Formatting results
form_comp_all <- formatBayesian(mortality, res_comp_all, data_comp_all, "Combined_all")
form_sev_all <- formatBayesian(mortality, res_sev_all, data_sev_all, "Severity_all", fixed = TRUE)
form_comp_us <- formatBayesian(mortality, res_comp_us, data_comp_us, "Combined_us")
form_sev_us <- formatBayesian(mortality, res_sev_us, data_sev_us, "Severity_us", fixed = TRUE)
form_comp_nonus <- formatBayesian(mortality, res_comp_nonus, data_comp_nonus, "Combined_nonus")
form_sev_nonus <- formatBayesian(mortality, res_sev_nonus, data_sev_nonus, "Severity_nonus", fixed = TRUE)
form_comp_post <- formatBayesian(mortality, res_comp_post, data_comp_post, "Combined_post")
form_sev_post <- formatBayesian(mortality, res_sev_post, data_sev_post, "Severity_post", fixed = TRUE)
form_comp_pre <- formatBayesian(mortality, res_comp_pre, data_comp_pre, "Combined_pre")
form_sev_pre <- formatBayesian(mortality, res_sev_pre, data_sev_pre, "Severity_pre", fixed = TRUE)
form_san <- formatBayesian(mortality, res_san, data_san, "Sanatorium/hospital")
form_nosan <- formatBayesian(mortality, res_nosan, data_nosan, "Non-Sanatorium")

## Saving results
save(form_comp_all, form_sev_all, form_comp_us, form_sev_us, form_comp_nonus, form_sev_nonus,
     form_comp_post, form_sev_post, form_comp_pre, form_sev_pre, form_san, form_nosan,
     file = "R/bayesian_mortality.RData")



#### Diagnostic Plots ------------------------------------------------------------------------------

## All Studies
png("Figures/xyplot_comp_all.png")
xyplot(eval_comp_all)
dev.off()
png("Figures/autocorr_comp_all.png")
autocorr.plot(eval_comp_all)
dev.off()

png("Figures/xyplot_sev_all.png")
xyplot(eval_sev_all)
dev.off()
png("Figures/autocorr_sev_all.png")
autocorr.plot(eval_sev_all)
dev.off()


## US Studies
png("Figures/xyplot_comp_us.png")
xyplot(eval_comp_us)
dev.off()
png("Figures/autocorr_comp_us.png")
autocorr.plot(eval_comp_us)
dev.off()

png("Figures/xyplot_sev_us.png")
xyplot(eval_sev_us)
dev.off()
png("Figures/autocorr_sev_us.png")
autocorr.plot(eval_sev_us)
dev.off()


## Non-US Studies
png("Figures/xyplot_comp_nonus.png")
xyplot(eval_comp_nonus)
dev.off()
png("Figures/autocorr_comp_nonus.png")
autocorr.plot(eval_comp_nonus)
dev.off()

png("Figures/xyplot_sev_nonus.png")
xyplot(eval_sev_nonus)
dev.off()
png("Figures/autocorr_sev_nonus.png")
autocorr.plot(eval_sev_nonus)
dev.off()


## US pre-1930s Studies
png("Figures/xyplot_comp_pre.png")
xyplot(eval_comp_pre)
dev.off()
png("Figures/autocorr_comp_pre.png")
autocorr.plot(eval_comp_pre)
dev.off()

png("Figures/xyplot_sev_pre.png")
xyplot(eval_sev_pre)
dev.off()
png("Figures/autocorr_sev_pre.png")
autocorr.plot(eval_sev_pre)
dev.off()


## US post-1930s Studies
png("Figures/xyplot_comp_post.png")
xyplot(eval_comp_post)
dev.off()
png("Figures/autocorr_comp_post.png")
autocorr.plot(eval_comp_post)
dev.off()

png("Figures/xyplot_sev_post.png")
xyplot(eval_sev_post)
dev.off()
png("Figures/autocorr_sev_post.png")
autocorr.plot(eval_sev_post)
dev.off()


## Sanatorium stratified
png("Figures/xyplot_san.png")
xyplot(eval_san)
dev.off()
png("Figures/autocorr_san.png")
autocorr.plot(eval_san)
dev.off()

png("Figures/xyplot_nosan.png")
xyplot(eval_nosan)
dev.off()
png("Figures/autocorr_nosan.png")
autocorr.plot(eval_nosan)
dev.off()

