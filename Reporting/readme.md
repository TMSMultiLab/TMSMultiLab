# Analysis of 32 meta-analyses which assessed reporting in TMS studies

This dataset and analysis is part of the TMS-RAT project. It forms the 'scoping review' of that work. Overall, the project is up-dating and significantly expanding the TMS reporting methods checklist of [Chipchase et al. 2012](https://doi.org/10.1016/j.clinph.2012.05.003).

The raw data files below were compiled from 32 studies which reported a meta-analysis of TMS studies, and reported 4 or more reporting checklist scores. A total of 681 studies were rated by the 32 individual reports. The individual study and rating data was extracted from the 32 studies by two independent raters (the first rater did all 32 studies), and discrepancies were resolved by the first rater.


## Raw data

<code>Chipchase_items_raw.csv</code> Raw data file containing the following columns:

- META - Meta-analysis or systematic review article (one of 32)
- PMID - PubMed ID of the included article (if available; if not: the DOI)
- AUTHOR - First author of the original, individual study
- YEAR - Year of publication, from PubMed data
- AgeR - The 'Age' item of the Chipchase questionnaire, Reported
- AgeC - The 'Age' item of the Chipchase questionnaire, Controlled
- ...UncC - the other 29 items of the Chipchase questionnaire, with Reported and Controlled for each
- TOT - Total score, in percent, across Reported and Controlled
- NUMR - Number of Items with Reported data
- NUMC - Number of Items with Controlled data
- SUMR - Sum of Items Reported
- SUMC - Sum of Items Controlled
- PERCR - Percentage of Items Reported
- PERCC - Percentage of Items Controlled
- R1R - Rater 1 rating of Reported
- R1C - Rater 1 rating of Controlled
- R2R - Rater 2 rating of Reported
- R2C - Rater 2 rating of Controlled
- PERCAG - Percentage agreement between raters
- PERCCH - Percentage chance agreement
- AC1 - Gwet's AC1 inter-rater reliability
- AC95U - Gwet's AC1 upper 95% confidence limit
- AC95L - Gwet's AC1 lower 95% confidence limit
- KAPPA - Kappa inter-rater reliability
- KAPU - Kappa upper 95% confidence limit
- KAPL - Kappa lower 95% confidence limit

<code>Chipchase_items_means.csv</code> - Raw data file containing the same data as the indidual Items above, but with only one row per meta-analysis or systematic review, and all data are the mean reporting rates across all included studies.

<code>Chipchase_studies.csv</code> - Raw data file containing the same data as the individual Items above, but with data in only one, two or three columns, giving the study mean Reported and/or Controlled score.


## Analysis

<code>Chipchase_reporting.m</code> - the analysis file loads the data, removes duplicates, downloads meta-data from PubMed, analyses and plots graphs to illustrate the data.


## Outputs

<code>Chipchase_pmids.mat</code> - Matlab data structure containing the bibliographic details of all 681 articles.

<code>Chipchase_item_reporting.png</code> - Mean and bootstrapped 95% CIs for Reported (in red) and Controlled Items across all studies with individual Item-level reporting (from Chipchase_items.csv).

<code>Chipchase_study_reporting.png</code> - Mean study-level Reported (red) and Controlled (blue) Items for all individual studies ordered by publication year (x-axis) in the top two panels, and for the mean (and 95% CI) of all studies grouped by year in the bottom two panels.

<code>Chipchase_study_reporting_journals.png</code> - Mean study-level Reported Items for individual studies ordered by publication year (x-axis) and grouped by colour for different journals that the articles were published in. Solid lines give the best-fit regression for Reported over Year, dotted lines give the 95% CI for the fitted line.

<code>Chipchase_study_reporting_journal_means.png</code> Mean and 95% CIs Reported Items over all years separately by journal (x-axis).
