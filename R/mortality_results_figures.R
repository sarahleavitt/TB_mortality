#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Mortality Analysis

##############################################################################
# This program creates figures of the mortality analysis results
##############################################################################

options(scipen=999)
options(digits = 10)

rm(list = ls())
source("R/utils.R")
reload_source()

#Reading in the study_id correspondence table
studyid <- read.csv("data/study_id.csv")

#Reading in analysis results
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

# Combined plot for all stratifications
surv_all <- bind_rows(form_pre$surv_dens, form_post$surv_dens,
                      form_namerica$surv_dens, form_europe$surv_dens,
                      form_yessan$surv_dens, form_nosan$surv_dens) %>%
  mutate(label = factor(label, levels = c("Pre-1930", "North America",
                                          "Sanatorium/hospital",
                                          "Post-1930", "Europe",
                                          "Non-Sanatorium")))

ind_all <- bind_rows(form_pre$ind_surv, form_post$ind_surv,
                     form_namerica$ind_surv, form_europe$ind_surv,
                     form_yessan$ind_surv, form_nosan$ind_surv) %>%
  mutate(label = factor(label, levels = c("Pre-1930", "North America",
                                          "Sanatorium/hospital",
                                          "Post-1930", "Europe",
                                          "Non-Sanatorium")))
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

ggsave("Figures/survival_curve_all.png", width = 6.5, height = 5)


##### Forest plots -------------------------------------------------------------

# Forest plot for combined model, sorted by date, stratifed by study location,
# colored by treatment location
pred_plot_all <- form_comb$pred_comb %>%
  left_join(studyid) %>%
  arrange(desc(enrollment_start)) %>%
  mutate(location = ifelse(is.na(location), "", location),
         location = factor(location, levels = c("North America",
                                                "Europe", "")),
         author_time = ifelse(first_author != "Overall",
                              paste0(first_author, ": ", time_period),
                              "Overall"),
         author_time = factor(author_time, levels = unique(author_time)))

ggplot(pred_plot_all %>% filter(value != "median"),
       aes(x = est, y = author_time, xmin = cilb, xmax = ciub)) +
  facet_grid(location ~ pred_label, scales = "free", space = "free") +
  geom_point(aes(color = sanatorium)) +
  geom_point(data = pred_plot_all %>% filter(shape == "Overall",
                                              value != "median"),
             color = 'black', shape = 18, size = 3) +
  geom_errorbar(aes(color = sanatorium), width = 0.5) +
  geom_errorbar(data = pred_plot_all %>% filter(shape == "Overall",
                                                 value != "median"),
                width = 0.5) +
  scale_x_continuous(name = "Survival Probability", limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  scale_color_discrete(name = "Treatment Location", na.translate = FALSE,
                       labels = c("Non-Sanatorium", "Sanatorium/hospital")) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom")

ggsave("Figures/forest_comb.png", width = 6.5, height = 4.5)



# Forest plot for stratifications each separate
pred_plot_time <- bind_rows(form_pre$pred_comb, form_post$pred_comb) %>%
  group_by(label) %>%
  arrange(desc(first_author)) %>%
  mutate(first_author = factor(first_author, levels = unique(first_author))) %>%
  mutate(labelf = ifelse(shape == "Overall" & label == "Pre-1930", " ",
                         ifelse(shape == "Overall" & label == "Post-1930", "  ",
                                label)),
         labelf = factor(labelf, levels = c("Pre-1930", " ", "Post-1930", "  ")))

pred_plot_loc <- bind_rows(form_namerica$pred_comb, form_europe$pred_comb) %>%
  group_by(label) %>%
  arrange(desc(first_author)) %>%
  mutate(first_author = factor(first_author, levels = unique(first_author))) %>%
  mutate(labelf = ifelse(shape == "Overall" & label == "North America", " ",
                         ifelse(shape == "Overall" & label == "Europe", "  ",
                                label)),
         labelf = factor(labelf, levels = c("North America", " ", "Europe", "  ")))

pred_plot_san <- bind_rows(form_yessan$pred_comb, form_nosan$pred_comb) %>%
  group_by(label) %>%
  arrange(desc(first_author)) %>%
  mutate(first_author = factor(first_author, levels = unique(first_author))) %>%
  mutate(labelf = ifelse(shape == "Overall" & label == "Sanatorium/hospital", " ",
                         ifelse(shape == "Overall" & label == "Non-Sanatorium", "  ",
                                label)),
         labelf = factor(labelf, levels = c("Sanatorium/hospital", " ",
                                            "Non-Sanatorium", "  ")))

plot_forest <- function(pred_plot){
  p <- ggplot(pred_plot %>% filter(value != "median"),
              aes(x = est, y = first_author, xmin = cilb, xmax = ciub)) +
  geom_point() +
  geom_point(data = pred_plot %>% filter(shape == "Overall",
                                             value != "median"),
             color = 'black', shape = 18, size = 3) +
  geom_errorbar(width = 0.5) +
  scale_x_continuous(name = "Survival Probability", limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  facet_grid(labelf ~ pred_label, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom")
  
  return(p)
}

p_time <- plot_forest(pred_plot_time)
p_loc <- plot_forest(pred_plot_loc)
p_san <- plot_forest(pred_plot_san)

ggsave("Figures/forest_time.png", p_time, width = 6.5, height = 4)
ggsave("Figures/forest_loc.png", p_loc, width = 6.5, height = 4)
ggsave("Figures/forest_san.png", p_san, width = 6.5, height = 4)

