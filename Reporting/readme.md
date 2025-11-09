# Analysis of 32 meta-analyses which assessed reporting in TMS studies

This dataset and analysis is part of the TMS-RAT project. It forms the 'scoping review' of that work. Overall, the project is up-dating and significantly expanding the TMS reporting methods checklist of [Chipchase et al. 2012](https://doi.org/10.1016/j.clinph.2012.05.003).

The raw data files below were compiled from 32 studies which reported a meta-analysis of TMS studies, and reported 4 or more reporting checklist scores. A total of 681 studies were rated by the 32 individual reports. The individual study and rating data was extracted from the 32 studies by two independent raters (the first rater did all 32 studies), and discrepancies were resolved by the first rater.


## Raw data

<code>Chipchase_items_raw.csv</code> This raw data file contains the following columns:

- META - Meta-analysis or systematic review article (one of 32)
- PMID - PubMed ID of the included article (if available; if not: the DOI)
- AUTHOR - First author of the original, individual study
- YEAR - Year of publication, from PubMed data
- AgeR - The 'Age' item of the Chipchase questionnaire, Reported
- AgeC - The 'Age' item of the Chipchase questionnaire, Controlled
- ...UncC - the other 29 items of the Chipchase questionnaire, with Reported and Controlled for each
- TOT
- NUMR
- NUMC
- SUMR
- SUMC
- PERCR
- PERCC
- R1R
- R1C
- R2R
- R2C
- PERCAG
- PERCCH
- AC1
- AC95U
- AC95L
- KAPPA
- KAPU
- KAPL

<code>Chipchase_items_means.csv</code>

<code>Chipchase_studies.csv</code>


## Analysis

<code>Chipchase_reporting.m</code> - the analysis file loads the data, removes duplicates, downloads meta-data from PubMed, analyses and plots graphs to illustrate the data.


## Outputs

<code>Chipchase_pmids.mat</code> - a matlab structure containing the bibliographic details of all 681 articles

<code>Chipchase_item_reporting.png</code>

<code>Chipchase_study_reporting.png</code>

<code>Chipchase_study_reporting_journals.png</code>

<code>Chipchase_study_reporting_journal_means.png</code>
