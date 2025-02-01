%% GET A PUBMED RECORD
function record=PubMedGet(pmid)
    if nargin~=1
        error('wrong number of arguments: only one numeric argument (PubMed ID) required');
    end
    if ~isnumeric(pmid)
        error('wrong data type: only one numeric argument (PubMed ID) required');
    else
        % get the web data
        url=['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=',int2str(pmid),'&retmode=xml'];
	pause(5); % wait 5 seconds between every query to avoid spamming PMID & getting kicked off the web
	raw=urlread(url);
	record=PubMedParse(raw);
    end
end