%% load chipchase reporting data and do analysis on it, to be published here: https://github.com/TMSMultiLab/TMSMultiLab/wiki/Reporting-checklist

% q 1-8     participants
% q 9-25    methods
% q 26-28   paired-pulse
% q 29-30   analysis

%% LOAD THE DATA___________________________________________________________
items_raw=readtable('Chipchase_items_raw.csv');                             % studies with individual items reported: per study and per item
studies=readtable('Chipchase_studies.csv');                                 % studies with overall reporting score only: per study
items_means=readtable('Chipchase_items_means.csv');                         % meta-studies with overall item score only: per item

%% REMOVE DUPLICATE STUDIES________________________________________________
items_raw_duplicates=[];
studies_duplicates=[];

% items_raw table
ix=isfinite(items_raw.PMID);                                                % which have PMIDS?
pmids=sortrows(items_raw.PMID(ix));                                         % list of PMIDS in order
duplicates=find(diff(pmids)==0);                                            % any duplicates?
for p=1:numel(duplicates)                                                   % for each duplicate PMID
    ix=find(items_raw.PMID==pmids(duplicates(p)));                          % which studies in table have same PMID?
    
    % save duplicates for later analysis
    items_raw_duplicates=[items_raw_duplicates;items_raw(ix,:)];            % append data to new table
    d=nanmean(table2array(items_raw(ix,5:end)));                            % get mean across valid inputs
    items_raw(ix(1),5:end) = array2table(d);                                % put the mean back into the first repeat
    items_raw(ix(2:end),:)=[];                                              % remove the duplicated data
end
% studies table
ix=isfinite(studies.PMID);                                                  % which have PMIDS?
pmids=sortrows(studies.PMID(ix));                                           % list of PMIDS in order
duplicates=find(diff(pmids)==0);                                            % any duplicates?
for p=1:numel(duplicates)                                                   % for each duplicate PMID
    ix=find(studies.PMID==pmids(duplicates(p)));                            % which studies in table have same PMID?
    
    % save duplicates for later analysis
    studies_duplicates=[studies_duplicates;studies(ix,:)];                  % append data to new table
    
    d=nanmean(table2array(studies(ix,5:end)));                              % get mean across valid inputs
    studies(ix(1),5:end) = array2table(d);                                  % put the mean back into the first repeat
    studies(ix(2:end),:)=[];                                                % remove the duplicated data
end

%% PARAMETERS______________________________________________________________
plotcols={'r','b'};
jitter=[-0.25,0.25];
items=30;
its=10000; % iterations for bootstrap analysis

%% ADD TOTALS______________________________________________________________
TotR=nansum(table2array(items_raw(:,5:2:items.*2+4)),2);
TotC=nansum(table2array(items_raw(:,6:2:items.*2+4)),2);
items_raw=addvars(items_raw,TotR,TotC);

%% ADD COUNTS______________________________________________________________
CouR=sum(isfinite(table2array(items_raw(:,5:2:items.*2+4))),2);
CouC=sum(isfinite(table2array(items_raw(:,6:2:items.*2+4))),2);
items_raw=addvars(items_raw,CouR,CouC);

%% ADD PERCENTS____________________________________________________________
PerR=TotR./CouR;
PerC=TotC./CouC;
items_raw=addvars(items_raw,PerR,PerC);

%% EXTRACT CITATION DATA FROM PUBMED & SAVE (WARNING: DO NOT SPAM PUBMED! - 4-8 seconds between each query - )
if exist('Chipchase_pmids.mat','file')==0                                   % if the PMID data have not yet been extracted...
    ix=isfinite(items_raw.PMID);                                            % which have PMIDS?
    pmids=sortrows(items_raw.PMID(ix));                                     % list of PMIDS from items_raw table, in order
    ix=isfinite(studies.PMID);                                              % which have PMIDS?
    pmids=sortrows([pmids;sortrows(studies.PMID(ix))]);                     % list of PMIDS from studies table, in order
    duplicates=diff(pmids)==0;                                              % which are duplicated
    pmids=pmids(~duplicates);                                               % remove the duplicates
    pubmeddata=struct();                                                    % for saving the pubmed data
    addpath('../Metaanalysis/');                                            % add path to TMSMultiLab PubMed functions
    
    % additional PMID to add
    pmids(size(pmids,1)+1)=37586677;

    for p=572%1:numel(pmids)                                                    % for each unique PMID
        disp([' Retrieving record ',int2str(p),'...']);                     % keep us informed
        raw=PubMedGet(pmids(p));                                            % retrieve the PubMed XML record
	pubmeddata(p)=PubMedParse(raw);                                     % parse the PubMed XML data into a structure
    end
    save('Chipchase_pmids.mat','pubmeddata');
else
    load('Chipchase_pmids.mat');                                            % load the previously-extracted data
end

%% PROPORTION OF STUDIES THAT MEET EACH CRITERION__________________________
figure(1);
hold on;
for q=1:items
    wiki=['Q',int2str(q),': '];                                             % build a text string for updating the wiki page
    for r=1:2
        % items from the items_raw table
        s=(q-1).*2+r+4;                                                     % column where the data is
        ix=isfinite(table2array(items_raw(:,s)));                           % rows with valid data
        n=sum(ix);                                                          % total studies assessed
        t=sum(table2array(items_raw(ix,s)));                                % total studies meeting the criterion

	% items from items_means table
	s=(q-1).*2+r+2;                                                     % column where the data is
	m=sum(isfinite(table2array(items_means(:,s))));                     % how many additional meta-analyses to add?
	n2=nansum(table2array(items_means(m,2)));                           % total N in these meta-analyses
	t2=nansum(round((table2array(items_means(m,2)).*table2array(items_means(m,s)))/100));% total studies meeting the criterion
	
        % bootstrap a 95% CI
        sample=[ones(round(t+t2),1);zeros(round(n+n2-t-t2),1)];             % put all the zeros and ones together
        bn=numel(sample);                                                   % total sample size
        bs=nan(its,1);                                                      % empty variable for the bootstrap
        for i=1:its                                                         % for each iteration
            bs(i)=mean(sample(ceil(rand(bn,1).*bn)));                       % get mean of a random sample with replacement
        end
        bs=sortrows(bs);                                                    % arrange from small to large
        plot(q+jitter(r),(t+t2)./(n+n2),[plotcols{r},'o'],'MarkerFaceColor',plotcols{r});% plot the mean data
        plot([q+jitter(r),q+jitter(r)],[bs(250),bs(9750)],[plotcols{r},'-'],'LineWidth',2);% plot the 95% bootstrapped confidence interval
	
	% build string for reporting on wiki page
	wiki=[wiki,' **',num2str(100.*(t+t2)./(n+n2),3),'**<br><br>{',num2str(100.*bs(250),3),',&nbsp;',num2str(100.*bs(9750),3),'}|'];
    end
    disp([wiki,'|']);
end
plot([8.5,8.5],[0,1],'-','Color',[0.6,0.6,0.6],'LineWidth',2);
plot([25.5,25.5],[0,1],'-','Color',[0.6,0.6,0.6],'LineWidth',2);
plot([28.5,28.5],[0,1],'-','Color',[0.6,0.6,0.6],'LineWidth',2);
text(0.75,0.95,'Participants','FontSize',16);
text(8.75,0.95,'Methods','FontSize',16);
text(25.75,0.95,{'Paired-','pulse'},'FontSize',16);
text(28.75,0.95,'Analysis','FontSize',16);
xlabel('Item number');
ylabel('Proportion meeting criterion');
title('Proportion of studies meeting Chipchase criteria');
axis([0,31,0,1]);
set(gcf,'Position',[1,50,1200,600]);
set(gca,'FontSize',20);
print('Chipchase_item_reporting.png','-dpng');
close(1);

%% TOTAL REPORTING AND CONTROLLING SCORES OVER TIME________________________
% summary data from items_raw table
d=[items_raw.PerR,items_raw.PerC,items_raw.YEAR,items_raw.PMID];           % all the available data

% add data from studies table
ix=~isfinite(studies.R);                                                   % if 'R' data missing, use TOTAL data
tmp=studies;                                                               % copy temporarily
tmp.R(ix)=tmp.TOTAL(ix);                                                   % copy some TOTAL data into R column
d=[d;tmp.R./100,tmp.C./100,tmp.YEAR,tmp.PMID];                             % concatenate all R and T data with items_raw summaries
clear ix tmp;                                                              % remove temporary variables

%% CORRECT THE YEAR WITH THE PUBLICATION DATA FROM PUBMED DATA_____________
journals=strings(size(d,1),1);                                             % to hold the journal abbreviations
for n=1:size(d,1)                                                          % for each unique article
    if isfinite(d(n,4))
        j=find([pubmeddata.PMID]==d(n,4));                                 % find in the pubmed table
        str=int2str(pubmeddata(j).Year);                                   % build a date string, from the year...
        if ~isempty(pubmeddata(j).Month)
            str=[str,'-',pubmeddata(j).Month];                             % and the month...
	    if ~isempty(pubmeddata(j).Day)
	        str=[str,'-',int2str(pubmeddata(j).Day)];                  % and the day
	    else
	        str=[str,'-15'];                                           % assume middle of month if missing day
	    end
        else
            str=[str,'-July-01'];                                          % assume middle of year if missing month & day
        end
        d(n,3)=datenum(datetime(str))./365.25;                             % convert to fractional year
        journals{n}=pubmeddata(j).Jrnl;                                    % list of journal abbreviations
    end
end

% fix some missing journals data
js=unique(journals);                                                       % list of unique journals
j=numel(js);
journals_unique=cell(j,4);                                                 % Journal name, N, mean, SD
% get N and mean reporting per journal
for j=1:size(journals_unique,1)
    journals_unique{j,1}=js(j,:);
    disp(journals_unique{j,1});
    ix=find(strcmp(journals,journals_unique{j,1})==1);
    journals_unique{j,2}=numel(ix);
    journals_unique{j,3}=nanmean(d(ix,1));                                 % mean reporting
    journals_unique{j,4}=nanstd(d(ix,1));                                  % SD reporting
end

% index to valid data
clear ix;
ix(:,1)=isfinite(d(:,1));
ix(:,2)=isfinite(d(:,2));
xpred=linspace(min(d(:,3)),max(d(:,3)))';                                  % 100 equally spaced points across years

% plot the data
figure(2);
for r=1:2
    % plot raw data and regression fit
    subplot(2,2,r);
    hold on;
    plot([2012,2012],[0,1],'-','Color',[0.6,0.6,0.6],'LineWidth',1);
    plot(d(ix(:,r),3),d(ix(:,r),r),[plotcols{r},'o'],'MarkerFaceColor',plotcols{r},'MarkerSize',3);

    % fit linear model
    mdl=fitlm(d(ix(:,r),3),d(ix(:,r),r));
    [ypred,ci]=predict(mdl,xpred);
    plot(xpred,ypred,[plotcols{r},'-'],'Linewidth',2);
    plot(xpred,ci(:,1),[plotcols{r},'--']);
    plot(xpred,ci(:,2),[plotcols{r},'--']);
    
    text(1986,1,['R^2=',num2str(mdl.Rsquared.Adjusted,3)],'Color',plotcols{r},'FontSize',16);
    axis([1985,2025,0,1]);
    set(gca,'FontSize',20);
    
    % plot means
    subplot(2,2,r+2);
    hold on;
    plot([2012,2012],[0,1],'-','Color',[0.6,0.6,0.6],'LineWidth',1);   
    for y=1985:2025
        yx=logical(ix(:,r).*d(:,3)>=y & ix(:,r).*d(:,3)<y+1);               % index to studies in this year
        N=sum(yx);
        if N>2
            M=nanmean(d(yx,r));                                             % mean across all studies for this year
            SD=nanstd(d(yx,r));                                             % SD across all studies for this year
            SE=SD./sqrt(N);
            CI=SE.*tinv(.975,N-1);
            plot(y,M,[plotcols{r},'s'],'MarkerSize',8);
            plot([y,y],[M-CI,M+CI],[plotcols{r},'-'],'LineWidth',1);
        end
    end
    axis([1985,2025,0,1]);
    set(gca,'FontSize',20);
end
subplot(2,2,1);
ylabel('Proportion criteria met');
title('Reported');
subplot(2,2,2);
title('Controlled');
subplot(2,2,3);
xlabel('Year');
ylabel('Proportion criteria met');
subplot(2,2,4);
xlabel('Year');
set(gcf,'Position',[1,50,1200,600]);
print('Chipchase_study_reporting.png','-dpng');
close(2);


%% REPEAT, THIS TIME FOR DIFFERENT JOURNALS

%% find a journal name
% find(strcmp(journals,"Clin Neurophysiol"))                               % in string array of same size as the data
% find(strcmp(cellstr([journals_unique{:,1}]),'Clin Neurophysiol')==1)	   % in cell array with means across journals list

% sort journals by most-represented first
jix=nan(size(journals_unique,1),2);                                        % empty array
jix(:,2)=[journals_unique{:,2}];                                           % number of articles per journal
jix(:,1)=1:size(journals_unique,1);                                        % index to journals_unique
jix=sortrows(jix,2,'descend');                                             % list of most-represented journals first
plotcols={'r','b','g','c','k'};
journals_unique{1,1}="Other";                                              % missing name is 2nd-most common
journals(find(journals==""))="Other";                                      % replace missing journals

figure(3);
hold on;
plot([2012,2012],[0,1],'-','Color',[0.6,0.6,0.6],'LineWidth',1);
n=0;
for j=[1:5]
    k=mod(j-1,5)+1;
    l=ceil(j./5);
    ix=find(strcmp(journals,journals_unique{jix(j),1})==1);                % list of articles from this journal
    plot(d(ix,3),d(ix,1),[plotcols{k},'o'],'MarkerFaceColor',plotcols{k},'MarkerSize',5);
    mdl=fitlm(d(ix,3),d(ix,1));
    xpred=linspace(min(d(ix,3)),max(d(ix,3)))';
    [ypred,ci]=predict(mdl,xpred);
    plot(xpred,ypred,[plotcols{k},'-'],'Linewidth',2);
    plot(xpred,ci(:,1),[plotcols{k},':'],'Linewidth',1);
    plot(xpred,ci(:,2),[plotcols{k},':'],'Linewidth',1);
    label=strjoin([strip(journals_unique{jix(j),1}),', N=',int2str(journals_unique{jix(j),2})],'');
    text(1986,1-(j./20),label,'color',plotcols{k},'FontSize',16);
end
axis([1985,2025,0,1]);
set(gca,'FontSize',20);
xlabel('Year');
ylabel('Proportion criteria met');
set(gcf,'Position',[1,50,1200,600]);
print('Chipchase_study_reporting_journals.png','-dpng');
close(3);

% REPEAT FOR MEANS OF TOP 25 JOURNALS
figure(4);
hold on;
n=0;
for j=[1:25]
    k=mod(j-1,5)+1;
    l=ceil(j./5);
    ix=find(strcmp(journals,journals_unique{jix(j),1})==1);                % list of articles from this journal
    M=mean(d(ix,1));
    S=std(d(ix,1));
    plot(j,M,[plotcols{k},'o'],'MarkerFaceColor',plotcols{k},'MarkerSize',5);
    plot([j,j],[M-S,M+S],'-','Color',plotcols{k});
end
xticks(1:25);
xticklabels([journals_unique{jix([1:25],1),1}]);
set(gca,'FontSize',12);
ylabel('Proportion criteria met');
set(gcf,'Position',[1,50,1200,600]);
print('Chipchase_study_reporting_journal_means.png','-dpng');
close(4);
