# Constructing an open-source test dataset for MEP analyses

This directory will contain a series of scripts to retrieve and reformat open datasets for inclusion into a large, multi-lab open dataset.

The dataset will be used to test TMS, EMG, and MEP analysis code.

## Guidance

We don't want to store the data itself. We want to provide scripts that will allow users to download the data (semi-automatically), and recreate the full test dataset on their local machine.

One script per dataset.

For the first dataset:

* URL to download the data from
* Description of the dataset (plain text, context)
* Details of methods and analysis
* Data file (2D array, with one EMG epoch per row and one column per sample, e.g., 10000 rows and 2000 columns)
* Metadata file (2D array, with one row per row of the data file, and one column per metadata variable e.g., 10000 rows and 20 columns)
* Subjects file (2D array, with one row per subject and one column per subject metadata variable)

Future plans:
* Convert all data to BIDS format
