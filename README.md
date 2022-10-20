# Analysis of Pre-chemotherapy TB Studies

This directory contains the code and data to run analyses and produce results for 
TITLE by Rodriguez, CA and Leavitt SV, et al.
This paper describes a meta-analysis of disease prognosis (including survival
and self-cure rates) of pre-chemotherapy tuberculosis studies.

## Data

### cure_data.csv

This is the individual-level self-cure data created from pre_chemo_data.xlsx by 
data_prep.R.


### cure_data_all.csv

LAURA REMOVE OR DESCRIBE


### cureDataSummary.csv

LAURA REMOVE OR DESCRIBE (if keeping, change name to cure_data_summary.csv to be consistent)


### mortality_data.csv

This is the individual-level mortality data created from pre_chemo_data.xlsx by 
data_prep.R.


### mortality_results_table.csv

This is a table of the results of the mortality analysis,


### NaturalRecovaryData.csv

LAURA REMOVE OR DESCRIBE (if keeping, change name to natural_recovery_data.csv to be consistent)


### pre_chemo_data.xlsx

This is an Excel spreadsheet with the cleaned, extracted life-table data for each of 
the studies. The tabs are labeled with the numeric ID given to each study.


### study_id.csv

This table details the concordance between the numeric study IDs and the papers they 
refer to (first author and year) as well as other information about each study,


### study_severity.csv

This is a table of the disease severity distribution for each study that has that information.


### summary_severity.csv

This is a table of the disease severity distribution for sanatorium/hospital studies
verses non-sanatorium studies.

***

## Scripts

### data_prep.R

This script takes the Excel spreadsheet with the extracted life-table data 
(pre_chemo_data.xlsx) and formats it into the individual-level data for both 
mortality (mortality_data.csv) and self-cure (cure_data.csv) to be used in 
analysis. LAURA ADD WHAT YOU CREATE HERE


### mortality_analysis.R

This script runs the Bayesian mortality TB-survival analysis for the complete
and stratified models and saves the output in an R workspace (bayesian_raw.RMD).


### mortality_results_format.R

This script takes the results of the analyses and formats and saves the 
output in an R workspace (bayesian_clean.RMD). It also saves the diagnostic plots.


### mortality_results_tables.R

This script creates the tables for the main text and supplement.


### mortality_results_figures.R

This script creates the figures for the main text and supplement.


### plot_all_cure_data.R

LAURA REMOVE OR DESCRIBE


### plotAllCureData.html

LAURA REMOVE OR DESCRIBE (if keeping, change name to plot_all_cure_data.html for consistency).


### plotAllCureData.Rmd

LAURA REMOVE OR DESCRIBE (if keeping, change name to plot_all_cure_data.Rmd for consistency).


### Ragonnet calculations.R

LAURA REMOVE OR DESCRIBE (if keeping, change name to ragonnet_calculations.R for consitency)


### severity_distribution.R

This script finds the disease severity distribution comparing sanatorium/hospital studies to
non-sanatorium studies


### utils.R

This script contains many functions called by the other programs to format the data
for analysis and results for presentation.





