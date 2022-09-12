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
load('R/bayesian_mortality.RData')



#### Summary Survival Curves -----------------------------------------------------------------------

plot_sum_curves <- function(complete_model, stratified_model){
  
  #Survival curve for complete model
  p1 <- ggplot(complete_model$surv_dens) +
    geom_line(aes(x = x, y = surv),
              color = "black", size = 1, linetype = "solid") +
    geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
                stat = "identity", linetype = 0, alpha = 0.25, na.rm = TRUE) +
    scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
    scale_x_continuous(name = "Years", limits = c(0, 30)) +
    theme_bw()
  
  
  #TB survival for fixed effect model
  p2 <- ggplot(stratified_model$surv_dens) +
    geom_line(aes(x = x, y = surv, color = severity),
              size = 1, linetype = "solid") +
    geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub, fill = severity),
                stat = "identity", linetype = 0, alpha = 0.15, na.rm = TRUE) +
    scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
    scale_x_continuous(name = "Years", limits = c(0, 30)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual("", values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                      "Far advanced" = "firebrick2", "Unknown" = "grey50")) +
    scale_fill_manual("",
                      values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                 "Far advanced" = "firebrick2", "Unknown" = "grey50"))
  
  p_comb <- arrangeGrob(p1, p2, nrow = 2)
  return(p_comb)
}

#All studies
p_sum_all <- plot_sum_curves(form_comp_all, form_sev_all)
ggsave("Figures/summary_curves_all.png", p_sum_all, width = 5.5, height = 8)

#US studies
p_sum_us <- plot_sum_curves(form_comp_us, form_sev_us)
ggsave("Figures/summary_curves_us.png", p_sum_us, width = 5.5, height = 8)

#Non-US studies
p_sum_nonus <- plot_sum_curves(form_comp_nonus, form_sev_nonus)
ggsave("Figures/summary_curves_nonus.png", p_sum_nonus, width = 5.5, height = 8)

#US post-1930s studies
p_sum_post <- plot_sum_curves(form_comp_post, form_sev_post)
ggsave("Figures/summary_curves_post.png", p_sum_post, width = 5.5, height = 8)

#US pre-1930s studies
p_sum_pre <- plot_sum_curves(form_comp_pre, form_sev_pre)
ggsave("Figures/summary_curves_pre.png", p_sum_pre, width = 5.5, height = 8)



#### Individual Survival Curves --------------------------------------------------------------------

plot_ind_curves <- function(complete_model, stratified_model){
  
  #Survival curves for complete model
  p1 <- ggplot(complete_model$ind_surv) +
    geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
              size = 0.7, alpha = 0.3) +
    geom_line(data = complete_model$surv_dens, aes(x = x, y = surv),
              color = "black", size = 1, linetype = "longdash") +
    scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
    scale_x_continuous(name = "Years", limits = c(0, 30)) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual("Disease Severity",
                       values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                  "Far advanced" = "firebrick2", "Unknown" = "grey50"))
  
  #Survival curves for stratified model
  p2 <- ggplot(stratified_model$ind_surv) +
    geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
              size = 0.7, alpha = 0.3) +
    geom_line(data = stratified_model$surv_dens,
              aes(x = x, y = surv, color = severity),
              linetype = "longdash", size = 1) +
    scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
    scale_x_continuous(name = "Years", limits = c(0, 30)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual("", drop = FALSE,
                       values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                  "Far advanced" = "firebrick2", "Unknown" = "grey50"))
  
  p_comb <- arrangeGrob(p1, p2, nrow = 2)
  return(p_comb)
}

#All studies
p_ind_all <- plot_ind_curves(form_comp_all, form_sev_all)
ggsave("Figures/individual_curves_all.png", p_ind_all, width = 5.5, height = 8)

#US studies
p_ind_us <- plot_ind_curves(form_comp_us, form_sev_us)
ggsave("Figures/individual_curves_us.png", p_ind_us, width = 5.5, height = 8)

#Non-US studies
p_ind_nonus <- plot_ind_curves(form_comp_nonus, form_sev_nonus)
ggsave("Figures/individual_curves_nonus.png", p_ind_nonus, width = 5.5, height = 8)

#US post-1930s studies
p_ind_post <- plot_ind_curves(form_comp_post, form_sev_post)
ggsave("Figures/individual_curves_post.png", p_ind_post, width = 5.5, height = 8)

#US pre-1930s studies
p_ind_pre <- plot_ind_curves(form_comp_pre, form_sev_pre)
ggsave("Figures/individual_curves_pre.png", p_ind_pre, width = 5.5, height = 8)




#### All studies survival curves by category -------------------------------------------------------

ggplot(form_sev_all$ind_surv) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = category),
            size = 0.7, alpha = 0.3) +
  geom_line(data = form_sev_all$surv_dens,
            aes(x = x, y = surv, group = severity),
            linetype = "longdash", size = 1, color = "grey50") +
  facet_wrap(~severity) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("Type of Study",
                     values = c("Non-US" = "mediumturquoise", "US post-1930" = "mediumvioletred",
                                "US pre-1930" = "royalblue3")) +
  ggsave("Figures/survival_curves_category.png", width = 8, height = 4.5)



##### Forest plots --------------------------------------------------------------------------------

pred_plot_all <- form_comp_all$pred_comb %>%
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



#### Sanatorium Sensitivity Analysis----------------------------------------------------------------

ggplot(bind_rows(form_san$ind_surv, form_nosan$ind_surv)) +
  facet_wrap(~label, nrow = 2) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
            size = 1, alpha = 0.3) +
  geom_line(data = bind_rows(form_san$surv_dens, form_nosan$surv_dens),
            aes(x = x, y = surv),
            color = "black", size = 1, linetype = "longdash") +
  geom_smooth(data = bind_rows(form_san$surv_dens, form_nosan$surv_dens),
              aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
              stat = "identity", linetype = 0, alpha = 0.4) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("", values = c("Minimal" = "seagreen", "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2", "Unknown" = "grey50")) +
  ggsave("Figures/sanatorium_curves.png", width = 5.5, height = 8)




