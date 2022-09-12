# Analysis of Pre-chemotherapy TB Studies

This directory contains the code and data to run analyses and produce results for 
TITLE by Rodriguez, CA and Leavitt SV, et al.
This paper describes a meta-analysis of disease prognosis (including survival
and self-cure rates) of pre-chemotherapy tuberculosis studies stratified by disease
severity.

## Data

### pre_chemo_data.xlsx

This is an Excel spreadsheet with the cleaned, extracted life-table data for each of 
the studies. The tabs are labeled with the numeric ID given to each study.

### study_id.csv

This table details the concordance between the numeric study IDs and the papers they 
refer to (first author and year) as well as the category of each study (US pre-1930s,
US post-1930s, Non-US).

### mortality_data.csv

This is the individual-level mortality data created from pre_chemo_data.xlsx by 
data_prep.R.

### cure_data.csv

This is the individual-level self-cure data created from pre_chemo_data.xlsx by 
data_prep.R.


***

## Scripts

### utils.R

This script contains many functions called by the other programs to format the data
for analysis and results for presentation.

### data_prep.R

This script takes the Excel spreadsheet with the extracted life-table data 
(pre_chemo_data.xlsx) and formats it into the individual-level data for both 
mortality (mortality_data.csv) and self-cure (cure_data.csv) which will be used in 
analysis.


### mortality_functions.R

This script contains the functions to run all of the Bayesian mortality survival
models: TB-specific mortality both with (stratified model) and without (complete model)
a fixed effect for disease severity. These two models was used on all of the
studies combined as well as stratified by study type (US pre-1930s, US post-1930s,
Non-US). The complete model was used on a sensitivity analysis of santorium/hospital
verses non-sanatorium studies.

### mortality_analysis_all.R

This script runs and saves the results of the Bayesian mortality TB-survival analysis
for the complete and stratified model for all of the studies.

### mortality_analysis_us.R

This script runs and saves the results of the Bayesian mortality TB-survival analysis
for the complete and stratified model for all US studies.

### mortality_analysis_nonus.R

This script runs and saves the results of the Bayesian mortality TB-survival analysis
for the complete and stratified model for the non-US studies.

### mortality_analysis_post.R

This script runs and saves the results of the Bayesian mortality TB-survival analysis
for the complete and stratified model for the US post-1930s studies.

### mortality_analysis_pre.R

This script runs and saves the results of the Bayesian mortality TB-survival analysis
for the complete and stratified model for the US pre-1930s studies.

### mortality_analysis_all.R

This script runs and saves the results of the Bayesian mortality TB-survival analysis
for the complete model stratified by whether the studies were from santorium/hospital
or not.


### mortality_results_format.R

This script takes the results of the all of the analyses and formats and saves the 
output in an R workspace (bayesian_mortality.RMD) 

### mortality_results_tables.R

This script creates the tables for the main text and supplement.

### mortality_results_figures.R

This script creates the figures for the main text and supplement.



### cure_analysis.R

This script runs the Bayesian logistic regression model for self-cure and extracts the data from the results for the table used in the manuscript.


LAURA - ADD INFO ABOUT WHAT YOU ADDED
