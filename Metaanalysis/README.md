# Meta-analysis code and data

## PubMedGet.m
Pass a pubmed ID (PMID) to the script and it retrieves the XML record from the PubMed server. Includes a random 4-8 second delay with every call to reduce demand on the pubmed server).

<code>raw=PubMedGet(pmid);</code>

## PubMedParse.m
Pass the results of <code>PubMedGet.m</code> in to get a record for one article.

<code>record=PubMedParse(raw);</code>
