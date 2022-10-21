#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Mortality Analysis

##############################################################################
# This program finds the disease severity distribution comparing 
# sanatorium/hospital studies to non-sanatorium studies 
##############################################################################

options(scipen=999)
options(digits = 10)
set.seed(150183)

rm(list = ls())
source("R/utils.R")
reload_source()

#Reading in individual mortality data for modeling
mortality <- read.csv("data/mortality_data.csv")

#Reading in the study_id correspondence table
studyid <- read.csv("data/study_id.csv")


#### Table of severity distribution by study

severity <- mortality %>% 
  filter(severity != "Unknown") %>%
  mutate(sanatorium = ifelse(sanatorium == "Yes", "Sanatorium/hospital",
                             "Non-sanatorium"),
         severity = factor(severity,
                           levels = c("Minimal", "Moderate", "Advanced")))

study_severity <- severity %>%
  group_by(study_id) %>%
  summarize(n = n(),
            Minimal_n = sum(severity == "Minimal"),
            Moderate_n = sum(severity == "Moderate"),
            Advanced_n = sum(severity == "Advanced"),
            sanatorium = first(sanatorium),
            first_author = first(first_author)) %>%
  mutate(Minimal_p = round(100 * Minimal_n / n, 1),
         Moderate_p = round(100 * Moderate_n / n, 1),
         Advanced_p = round(100 * Advanced_n / n, 1),
         Minimal = paste0(Minimal_n, " (", Minimal_p, "%)"),
         Moderate = paste0(Moderate_n, " (", Moderate_p, "%)"),
         Advanced = paste0(Advanced_n, " (", Advanced_p, "%)")) %>%
  select(first_author, sanatorium, Minimal, Moderate, Advanced) %>%
  arrange(sanatorium)

write.csv(study_severity, "data/study_severity.csv")


#### Summary table of severity distribution by treatment location

summary_severity1 <- severity %>%
  group_by(sanatorium) %>%
  summarize(n = n(),
            Minimal_n = sum(severity == "Minimal"),
            Moderate_n = sum(severity == "Moderate"),
            Advanced_n = sum(severity == "Advanced")) %>%
  mutate(Minimal_p = round(100 * Minimal_n / n, 1),
         Moderate_p = round(100 * Moderate_n / n, 1),
         Advanced_p = round(100 * Advanced_n / n, 1))

summary_severity2 <- summary_severity1 %>%
  mutate(Minimal = paste0(Minimal_n, " (", Minimal_p, "%)"),
         Moderate = paste0(Moderate_n, " (", Moderate_p, "%)"),
         Advanced = paste0(Advanced_n, " (", Advanced_p, "%)")) %>%
  select(sanatorium, Minimal, Moderate, Advanced)

write.csv(summary_severity2, "data/summary_severity.csv")

chisq.test(table(severity$severity, severity$sanatorium))


#### Barplot of severity distribution by treatment location

severity_plot <- summary_severity1 %>%
  select(sanatorium, Minimal_p, Moderate_p, Advanced_p) %>%
  pivot_longer(!sanatorium, names_to = "severity", values_to = "percent") %>%
  mutate(severity = factor(severity,
                           levels = c("Minimal_p", "Moderate_p", "Advanced_p")))

ggplot(data = severity_plot, aes(x = sanatorium, y = percent,
                                 fill = severity)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ylab("Percentage of Patients") +
  xlab("Setting") +
  scale_fill_manual(name = "Disease Severity",
                      breaks = c("Minimal_p", "Moderate_p", "Advanced_p"),
                      labels = c("Minimal", "Moderately\nadvanced",
                                 "Far advanced"),
                      values = c("grey70", "grey40", "grey20")) +
  theme_bw()
  
ggsave("Figures/severity_barchart.png", dpi = 300, width = 5.5, height = 4)
          
