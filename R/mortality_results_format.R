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
load("R/bayesian_raw.RData")

## Formatting results
form_comb <- formatBayesian(mortality, res_comb, data_comb, "Combined")
form_pre <- formatBayesian(mortality, res_pre, data_pre, "Pre-1930")
form_post <- formatBayesian(mortality, res_post, data_post, "Post-1930")
form_namerica <- formatBayesian(mortality, res_namerica, data_namerica, "North America")
form_europe <- formatBayesian(mortality, res_europe, data_europe, "Europe")
form_yessan <- formatBayesian(mortality, res_yessan, data_yessan, "Sanatorium/hospital")
form_nosan <- formatBayesian(mortality, res_nosan, data_nosan, "Non-Sanatorium")

## Saving results
save(form_comb, form_pre, form_post, form_namerica, form_europe, form_yessan,
     form_nosan, file = "R/bayesian_clean.RData")



#### Diagnostic Plots ------------------------------------------------------------------------------

#Combined analysis
png("Figures/xyplot_comb.png")
xyplot(eval_comb)
dev.off()
png("Figures/autocorr_comb.png")
autocorr.plot(eval_comb)
dev.off()

#Time period stratified
png("Figures/xyplot_pre.png")
xyplot(eval_pre)
dev.off()
png("Figures/autocorr_pre.png")
autocorr.plot(eval_pre)
dev.off()
png("Figures/xyplot_post.png")
xyplot(eval_post)
dev.off()
png("Figures/autocorr_post.png")
autocorr.plot(eval_post)
dev.off()

#Location stratified
png("Figures/xyplot_namerica.png")
xyplot(eval_namerica)
dev.off()
png("Figures/autocorr_namerica.png")
autocorr.plot(eval_namerica)
dev.off()
png("Figures/xyplot_europe.png")
xyplot(eval_europe)
dev.off()
png("Figures/autocorr_europe.png")
autocorr.plot(eval_europe)
dev.off()

#Sanatorium stratified
png("Figures/xyplot_yessan.png")
xyplot(eval_yessan)
dev.off()
png("Figures/autocorr_yessan.png")
autocorr.plot(eval_yessan)
dev.off()
png("Figures/xyplot_nosan.png")
xyplot(eval_nosan)
dev.off()
png("Figures/autocorr_nosan.png")
autocorr.plot(eval_nosan)
dev.off()



