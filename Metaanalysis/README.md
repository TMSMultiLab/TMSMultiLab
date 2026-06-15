# Meta-analysis code and data

## PubMedGet.m
Pass a pubmed ID (pmid) to the script and it retrieves the XML record from the PubMed server. Includes a random 4-8 second delay with every call to reduce demand on the pubmed server).

`raw = PubMedGet(pmid);`

## PubMedParse.m
Pass the results of `PubMedGet.m` in to get a record for one article.

`record = PubMedParse(raw);`

## forest_plot.m
Create a forest plot for a meta analysis. Calls the `meta_analysis.m` function

`[options, stats] = forest_plot(Ms, SEs [, labels] [, options]);`

`Ms` - a list of study means

`SEs` - a list of study SEs (must be the same number of studies as `Ms`

`options` - many!

## meta_analysis.m
Run a meta-analysis on a set of study Means and SEs.

`stats = meta_analysis(Ms, SEs, [options]);`

`Ms` - a list of study means

`SEs` - a list of study SEs (must be the same number of studies as `Ms`

`options` - many!

## trim_and_fill.m
Run the trim and fill algorithm on a set of study means.


