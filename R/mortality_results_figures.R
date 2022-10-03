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

#Reading in individual mortality data and analysis results
mortality <- read.csv("data/mortality_data.csv")
load('R/bayesian_clean.RData')


#### Survival Curves: Combined -------------------------------------------------

ggplot(form_comb$surv_dens) +
  geom_line(aes(x = x, y = surv), color = "black", size = 1,
            linetype = "longdash") +
  geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
              stat = "identity", linetype = 0, alpha = 0.25, na.rm = TRUE) +
  geom_line(data = form_comb$ind_surv, aes(x = x, y = surv, group = study_id),
            size = 0.7, alpha = 0.2) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw()

ggsave("Figures/survival_curve_all.png", width = 6, height = 5.5)


#### Survival Curves: Stratified -----------------------------------------------

plot_surv_strata <- function(data1, data2){

  p <- ggplot(bind_rows(data1$surv_dens, data2$surv_dens)) +
    geom_line(aes(x = x, y = surv),
              color = "black", size = 1, linetype = "longdash") +
    geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
                stat = "identity", linetype = 0, alpha = 0.4) +
    facet_wrap(~label, nrow = 2) +
    geom_line(data = bind_rows(data1$ind_surv, data2$ind_surv),
              aes(x = x, y = surv, group = study_id),
              size = 0.6, alpha = 0.2) +
    scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
    scale_x_continuous(name = "Years", limits = c(0, 30)) +
    theme_bw()
}

surv_time <- plot_surv_strata(form_pre, form_post)
ggsave("Figures/survival_curve_time.png", surv_time, width = 4.5, height = 6)

surv_loc <- plot_surv_strata(form_namerica, form_europe)
ggsave("Figures/survival_curve_location.png", surv_loc, width = 5.5, height = 8)

surv_san <- plot_surv_strata(form_yessan, form_nosan)
ggsave("Figures/survival_curve_sanatorium.png", surv_san, width = 5.5, height = 8)




##### Forest plots --------------------------------------------------------------------------------

pred_plot_all <- form_comb_all$pred_comb %>%
  mutate(category = ifelse(is.na(category), "Overall", category),
         category = factor(category, levels = c("US post-1930", "US pre-1930",
                                                "Non-US", "Overall"))) %>%
  arrange(desc(category), desc(first_author)) %>%
  mutate(first_author = factor(first_author, levels=unique(first_author)))

pred_plot_sev <- form_sev_all$pred_comb %>%
  mutate(category = ifelse(is.na(category), "Overall", category),
         category = factor(category, levels = c("US post-1930", "US pre-1930", 
                                                "Non-US", "Overall"))) %>%
  arrange(desc(category), desc(first_author)) %>%
  mutate(first_author = factor(first_author, levels=unique(first_author)))


#TB survival for full model
ggplot(pred_plot_all %>% filter(value != "median"),
             aes(x = est, y = first_author, xmin = cilb, xmax = ciub,
                 shape = shape, color = category)) +
  facet_grid(severity ~ pred_label, scales = "free", space = "free") +
  geom_point(aes(color = category)) +
  geom_point(data = pred_plot_all %>% filter(shape == "Overall" & value != "median"),
             color = 'black', shape = 18, size = 3) +
  geom_errorbar(width = 0.5) +
  scale_x_continuous(name = "1-Year Survival Probability", limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom") +
  scale_color_manual("Type of Study",
                     values = c("Non-US" = "mediumturquoise", "US post-1930" = "mediumvioletred",
                                "US pre-1930" = "royalblue3", "Overall" = "black")) +
  scale_shape_discrete(guide = "none") +
  ggsave("Figures/forest_full.png", width = 7, height = 6)

#TB survival for stratified model
ggplot(pred_plot_sev %>% filter(value != "median"),
              aes(x = est, y = first_author, xmin = cilb, xmax = ciub, shape = shape,
                  color = category)) +
  geom_point(aes(color = category)) +
  geom_point(data = pred_plot_sev %>% filter(shape == "Overall" & value != "median"),
             color = 'black', shape = 18, size = 3) +
  geom_errorbar(width = 0.5) +
  scale_x_continuous(name = "1-Year Survival Probability", limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  facet_grid(severity ~ pred_label, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom") +
  scale_color_manual("Type of Study",
                     values = c("Non-US" = "mediumturquoise", "US post-1930" = "mediumvioletred",
                                "US pre-1930" = "royalblue3", "Overall" = "black")) +
  scale_shape_discrete(guide = "none") +
  ggsave("Figures/forest_stratified.png", width = 7, height = 6)








