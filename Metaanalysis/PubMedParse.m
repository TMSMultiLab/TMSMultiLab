%% PARSE A SINGLE PUBMED RECORD FROM RAW XML DATA
function record=PubMedParse(raw)
	if nargin~=1
		error('wrong number of arguments: one xml record is needed');
	end
	record=struct('PMID','','Journal','','Jrnl','','Volume','','Issue','','Year','','Month','','Day','','Title','','StartPage','','EndPage','','DOI','','Abstract','','Authors','','Language','');
    
	% PUMBED ID
	start=strfind(raw,'<PMID Version="1">');
	finish=strfind(raw,'</PMID>');
	record.PMID=str2num(raw(start+18:finish-1));

	% DOI
	start=strfind(raw,'<ArticleId IdType="doi">');
	if ~isempty(start)
		finish=strfind(raw,'</ArticleId>');
		record.DOI=raw(start(1)+24:finish(finish>start(1))-1);
	end

	% Journal string
	start=strfind(raw,'<Journal>');
	finish=strfind(raw,'</Journal>');
	str=raw(start+9:finish-1);

		% Journal long title
		start=strfind(str,'<Title>');
		finish=strfind(str,'</Title>');
		record.Journal=str(start+7:finish-1);
		
		% Journal short title
		start=strfind(str,'<ISOAbbreviation>');
		finish=strfind(str,'</ISOAbbreviation>');
		record.Jrnl=str(start+17:finish-1);
		
		% Volume
		start=strfind(str,'<Volume>');
		finish=strfind(str,'</Volume>');
		record.Volume=str(start+8:finish-1);
		
		% Issue
		start=strfind(str,'<Issue>');
		if ~isempty(start)
			finish=strfind(str,'</Issue>');
			record.Issue=str(start+7:finish-1);
		end
		
		% PubDate string
		start=strfind(str,'<PubDate>');
		finish=strfind(str,'</PubDate>');
		str2=str(start+9:finish-1);
    
			% Year
			start=strfind(str2,'<Year>');
			finish=strfind(str2,'</Year>');
			record.Year=str2num(str2(start+6:finish-1));

			% Month
			start=strfind(str2,'<Month>');
			if ~isempty(start)
				finish=strfind(str2,'</Month>');
				record.Month=str2(start+7:finish-1);
			end
			
			% Day
			start=strfind(str2,'<Day>');
			if ~isempty(start)
				finish=strfind(str2,'</Day>');
				record.Day=str2num(str2(start+5:finish-1));
			end

	% TITLE
	start=strfind(raw,'<ArticleTitle>');
	finish=strfind(raw,'</ArticleTitle>');
	record.Title=raw(start+14:finish-1);

	% Pagination string
	start=strfind(raw,'<Pagination>');
	finish=strfind(raw,'</Pagination>');
	str=raw(start+12:finish-1);

		% StartPage
		start=strfind(str,'<StartPage>');
		if ~isempty(start)
			finish=strfind(str,'</StartPage>');
			record.StartPage=str(start+11:finish-1);
		end

		% EndPage
		start=strfind(str,'<EndPage>');
		if ~isempty(start)
			finish=strfind(str,'</EndPage>');
			record.EndPage=str(start+9:finish-1);
		end
	
	% ABSTRACT
	start=strfind(raw,'<AbstractText>');
	if ~isempty(start)
		finish=strfind(raw,'</AbstractText>');
		record.Abstract=raw(start+14:finish-1);
	end
	
	% AUTHORS
	start=strfind(raw,'<AuthorList CompleteYN="Y">');
	finish=strfind(raw,'</AuthorList>');
	str=raw(start+27:finish-1);
	
		% N AUTHORS
		a=strfind(str,'<LastName>');
		record.Authors.N=numel(a);
		
		% LAST NAMES
		b=strfind(str,'</LastName>');
		for n=1:record.Authors.N
			record.Authors.LastName{n}=str(a(n)+10:b(n)-1);
		end

		% FORE NAMES
		a=strfind(str,'<ForeName>');
		b=strfind(str,'</ForeName>');
		for n=1:record.Authors.N
			record.Authors.ForeName{n}=str(a(n)+10:b(n)-1);
		end

		% INITIALS
		a=strfind(str,'<Initials>');
		b=strfind(str,'</Initials>');
		for n=1:record.Authors.N
			record.Authors.Initials{n}=str(a(n)+10:b(n)-1);
		end
		
		% AFFILILATIONS
		a=strfind(str,'<Affiliation>');
		b=strfind(str,'</Affiliation>');
		for n=1:numel(a)
			record.Authors.Affiliation{n}=str(a(n)+13:b(n)-1);
		end
		
	% LANGUAGE
	start=strfind(raw,'<Language>');
	finish=strfind(raw,'</Language>');
	record.Language=raw(start+10:finish-1);
		
end