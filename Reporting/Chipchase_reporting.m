%% load chipchase reporting data and do analysis on it, to be published here: https://github.com/TMSMultiLab/TMSMultiLab/wiki/Reporting-checklist

% q 1-8     participants
% q 9-25    methods
% q 26-28   paired-pulse
% q 29-30   analysis

%% LOAD THE ITEMS DATA_____________________________________________________
items_raw = readtable('Chipchase_items_raw.csv');                           % studies with individual items reported: per study and per item

%% REMOVE DUPLICATE STUDIES FROM ITEMS RAW TABLE___________________________
items_raw_duplicates = [];
ix = isfinite(items_raw.PMID);                                              % which have PMIDS?
pmids = sortrows(items_raw.PMID(ix));                                       % list of PMIDS in order
duplicates = find(diff(pmids) == 0);                                        % any duplicates? (N=7 at last check; 14 datasets copied to duplicates table)
for p = 1:numel(duplicates)                                                 % for each duplicate PMID
    ix = find(items_raw.PMID == pmids(duplicates(p)));                      % which studies in table have same PMID?
    
    % save duplicates for later analysis
    items_raw_duplicates = [items_raw_duplicates;items_raw(ix,:)];          % append data to new table
    d = nanmean(table2array(items_raw(ix,5:end)));                          % get mean across valid inputs
    items_raw(ix(1),5:end) = array2table(d);                                % put the mean back into the first repeat
    items_raw(ix(2:end),:) = [];                                            % remove the duplicated data
end

%% LOAD THE STUDIES DATA___________________________________________________
studies = readtable('Chipchase_studies.csv');                               % studies with overall reporting score only: per study

%% REMOVE DUPLICATE STUDIES FROM STUDIES TABLE_____________________________
studies_duplicates = [];
ix = isfinite(studies.PMID);                                                % which have PMIDS?
pmids = sortrows(studies.PMID(ix));                                         % list of PMIDS in order
duplicates = find(diff(pmids) == 0);                                        % any duplicates? (N=60 at last check; 120 datasets copied to duplicates table)
for p = 1:numel(duplicates)                                                 % for each duplicate PMID
    ix = find(studies.PMID == pmids(duplicates(p)));                        % which studies in table have same PMID?
    
    % save duplicates for later analysis
    studies_duplicates = [studies_duplicates;studies(ix,:)];                % append data to new table
    d = nanmean(table2array(studies(ix,5:end)));                            % get mean across valid inputs
    studies(ix(1),5:end) = array2table(d);                                  % put the mean back into the first repeat
    studies(ix(2:end),:) = [];                                              % remove the duplicated data
end

%% ANALYSIS & PLOTTING PARAMETERS__________________________________________
plotcols = {'r','b'};                                                       % reporting and controlled colours
jitter = [-0.125,0.125];                                                    % offset either side of item number
items = 30;                                                                 % total N items
its = 10000;                                                                % iterations for bootstrap analyses

%% ADD TOTAL REPORTED ITEMS________________________________________________
TotR = nansum(table2array(items_raw(:,5:2:items.*2+4)),2);                  % add up all 1s from the 5th column (1st item) to the 64th (30th item)
TotC = nansum(table2array(items_raw(:,6:2:items.*2+4)),2);                  % add up all 1s from the 6th column (1st item) to the 65th (30th item)
d = find(items_raw.SUMR-TotR~=0);                                           % checksum (should be the same, apart from items averaged across duplicates)
items_raw=addvars(items_raw,TotR,TotC);

%% ADD TOTAL ITEMS INCLUDED________________________________________________
CouR = sum(isfinite(table2array(items_raw(:,5:2:items.*2+4))),2);           % count all 0s and 1s from the 5th column (1st item) to the 64th (30th item)
CouC = sum(isfinite(table2array(items_raw(:,6:2:items.*2+4))),2);           % count all 0s and 1s from the 6th column (1st item) to the 65th (30th item)
items_raw = addvars(items_raw,CouR,CouC);

%% ADD PERCENT REPORTING___________________________________________________
PerR = TotR./CouR;                                                          % Total Reported divided by Total Assessed
PerC = TotC./CouC;                                                          % Total Controlled divided by Total Assessed
items_raw = addvars(items_raw,PerR,PerC);

%% EXTRACT CITATION DATA FROM PUBMED & SAVE (WARNING: DO NOT SPAM PUBMED! - 4-8 seconds between each query - )
if exist('Chipchase_pmids.mat','file') == 0                                 % if the PMID data have not yet been extracted...
    ix = isfinite(items_raw.PMID);                                          % which have PMIDS?
    pmids = sortrows(items_raw.PMID(ix));                                   % list of PMIDS from items_raw table, in order
    ix = isfinite(studies.PMID);                                            % which have PMIDS?
    pmids = sortrows([pmids;sortrows(studies.PMID(ix))]);                   % list of PMIDS from studies table, in order
    duplicates = diff(pmids)==0;                                            % which are duplicated
    pmids = pmids(~duplicates);                                             % remove the duplicates
    pubmeddata = struct();                                                  % for saving the pubmed data
    addpath('../Metaanalysis/');                                            % add path to TMSMultiLab PubMed functions

    for p = 1:numel(pmids)                                                  % for each unique PMID
        disp([' Retrieving record ',int2str(p),'...']);                     % keep us informed
        raw = PubMedGet(pmids(p));                                          % retrieve the PubMed XML record
	if p == 1
	    pubmeddata = PubMedParse(raw);                                  % parse the PubMed XML data into a structure
	else
	    pubmeddata(p) = PubMedParse(raw);                               % parse the PubMed XML data into a structure
	end
    end
    save('Chipchase_pmids.mat','pubmeddata');
else
    load('Chipchase_pmids.mat');                                            % load the previously-extracted data
end

%% LOAD THE ITEM MEANS DATA________________________________________________
items_means = readtable('Chipchase_items_means.csv');                       % meta-studies with overall item score only: per item

%% PROPORTION OF STUDIES THAT MEET EACH CRITERION__________________________
figure(1);
hold on;
for q = 1:items
    wiki = ['Q',int2str(q),': '];                                           % build a text string for updating the wiki page
    for r=1:2
        % items from the items_raw table
        s = (q-1).*2+r+4;                                                   % column where the data is
        ix = isfinite(table2array(items_raw(:,s)));                         % rows with valid data
        n = sum(ix);                                                        % total studies assessed
        t = sum(table2array(items_raw(ix,s)));                              % total studies meeting the criterion

	% items from items_means table
	s = (q-1).*2+r+4;                                                   % column where the data is
	m = isfinite(table2array(items_means(:,s)));                        % which additional meta-analyses to add?
	n2 = nansum(table2array(items_means(m,4)));                         % total N in these meta-analyses
	t2 = nansum(round((table2array(items_means(m,4)).*table2array(items_means(m,s)))/100));% total studies meeting the criterion
	
        % bootstrap a 95% CI
        sample = [ones(round(t+t2),1);zeros(round(n+n2-t-t2),1)];           % put all the zeros and ones together
        bn = numel(sample);                                                 % total sample size
        bs = nan(its,1);                                                    % empty variable for the bootstrap
        for i = 1:its                                                       % for each iteration
            bs(i) = mean(sample(ceil(rand(bn,1).*bn)));                     % get mean of a random sample with replacement
        end
        bs = sortrows(bs);                                                  % arrange from small to large
        plot(q+jitter(r),(t+t2)./(n+n2),[plotcols{r},'o'],'MarkerFaceColor',plotcols{r});% plot the mean data
        plot([q+jitter(r),q+jitter(r)],[bs(250),bs(9750)],[plotcols{r},'-'],'LineWidth',2);% plot the 95% bootstrapped confidence interval
	
	% build string for reporting on wiki page
	wiki = [wiki,' **',num2str(100.*(t+t2)./(n+n2),3),'**<br><br>{',num2str(100.*bs(250),3),',&nbsp;',num2str(100.*bs(9750),3),'}|'];
    end
    disp([wiki,'|']);
end
plot([8.5,8.5],[0,1],'-','Color',[0.6,0.6,0.6],'LineWidth',2);              % line between participant items and method items
plot([25.5,25.5],[0,1],'-','Color',[0.6,0.6,0.6],'LineWidth',2);            % line between method items and paired-pulse items
plot([28.5,28.5],[0,1],'-','Color',[0.6,0.6,0.6],'LineWidth',2);            % line between paired-pulse items and analysis items
text(0.75,1.05,'Participants','FontSize',16);                               % labels for the different sections
text(8.75,1.05,'Methods','FontSize',16);
text(25.75,1.05,{'Paired-','pulse'},'FontSize',16);
text(28.75,1.05,'Analysis','FontSize',16);
xticks([1:30]);
xticklabels({'','','','','5','','','','','10','','','','','15','','','','','20','','','','','25','','','','','30'});
xtickangle(0);
xlabel('Item number');
ylabel('Proportion meeting criterion');
%title('Proportion of studies meeting Chipchase criteria');
axis([0,31,0,1]);
set(gcf,'Position',[1,50,1200,600]);
set(gca,'FontSize',20);
print('Chipchase_item_reporting.png','-dpng');
close(1);

%% TOTAL REPORTING AND CONTROLLING SCORES OVER TIME_______________________
% summary data from items_raw table
d = [items_raw.PerR,items_raw.PerC,items_raw.YEAR,items_raw.PMID];         % all the available data

% add data from studies table
ix =~ isfinite(studies.PERCR);                                             % if 'R' data missing, use TOTAL data
tmp = studies;                                                             % copy temporarily
tmp.PERCR(ix) = tmp.TOT(ix);                                               % copy some TOTAL data into R column
d = [d;tmp.PERCR./100,tmp.PERCC./100,tmp.YEAR,tmp.PMID];                   % concatenate all R and T data with items_raw summaries
clear ix tmp;                                                              % remove temporary variables

%% CORRECT THE YEAR WITH THE PUBLICATION DATE FROM PUBMED DATA____________
journals = strings(size(d,1),1);                                           % to hold the journal abbreviations
for n = 1:size(d,1)                                                        % for each unique article
    if isfinite(d(n,4))
        j = find([pubmeddata.PMID] == d(n,4));                             % find in the pubmed table
        str = int2str(pubmeddata(j).Year);                                 % build a date string, from the year...
        if ~isempty(pubmeddata(j).Month)
            str = [str,'-',pubmeddata(j).Month];                           % and the month...
	    if ~isempty(pubmeddata(j).Day)
	        str = [str,'-',int2str(pubmeddata(j).Day)];                % and the day
	    else
	        str = [str,'-15'];                                         % assume middle of month if missing day
	    end
        else
            str = [str,'-July-01'];                                        % assume middle of year if missing month & day
        end
        d(n,3) = datenum(datetime(str))./365.25;                           % convert to fractional year
        journals{n} = pubmeddata(j).Jrnl;                                  % list of journal abbreviations
    end
end

%% FIX SOME MISSING JOURNALS DATA_________________________________________
js = unique(journals);                                                     % list of unique journals
j = numel(js);                                                             % how many are there?
journals_unique = cell(j,4);                                               % Journal name, N, mean, SD
% get N and mean reporting per journal
for j = 1:size(journals_unique,1)
    journals_unique{j,1} = js(j,:);
    disp(journals_unique{j,1});
    ix = find(strcmp(journals,journals_unique{j,1}) == 1);
    journals_unique{j,2} = numel(ix);
    journals_unique{j,3} = nanmean(d(ix,1));                               % mean reporting
    journals_unique{j,4} = nanstd(d(ix,1));                                % SD reporting
end

%% CREATE INDEX TO VALID DATA_____________________________________________
clear ix;
ix(:,1) = isfinite(d(:,1));                                                % valid reporting data
ix(:,2) = isfinite(d(:,2));                                                % valid controlled data
xpred = linspace(min(d(:,3)),max(d(:,3)))';                                % 100 equally spaced points across min to max years

%% PLOT THE DATA__________________________________________________________
figure(2);
for r = 1:2
    % plot raw data and regression fit
    subplot(2,2,r);
    hold on;
    plot([2012,2012],[0,1],'-','Color',[0.6,0.6,0.6],'LineWidth',1);
    plot(d(ix(:,r),3),d(ix(:,r),r),[plotcols{r},'o'],'MarkerFaceColor',plotcols{r},'MarkerSize',3);

    % fit linear model
    mdl = fitlm(d(ix(:,r),3),d(ix(:,r),r));
    [ypred,ci] = predict(mdl,xpred);
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
    for y = 1985:2025
        yx = logical(ix(:,r).*d(:,3)>=y & ix(:,r).*d(:,3)<y+1);             % index to studies in this year
        N = sum(yx);
        if N>2
            M = nanmean(d(yx,r));                                           % mean across all studies for this year
            SD = nanstd(d(yx,r));                                           % SD across all studies for this year
            SE = SD./sqrt(N);
            CI = SE.*tinv(.975,N-1);
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


%% REPEAT, THIS TIME FOR DIFFERENT JOURNALS_______________________________
%% find a journal name
% find(strcmp(journals,"Clin Neurophysiol"))                               % in string array of same size as the data
% find(strcmp(cellstr([journals_unique{:,1}]),'Clin Neurophysiol')==1)	   % in cell array with means across journals list

%% sort journals by most-represented first________________________________
jix = nan(size(journals_unique,1),2);                                      % empty array
jix(:,2) = [journals_unique{:,2}];                                         % number of articles per journal
jix(:,1) = 1:size(journals_unique,1);                                      % index to journals_unique
jix = sortrows(jix,2,'descend');                                           % list of most-represented journals first
plotcols = {'r','b','g','c','k'};
journals_unique{1,1} = "Other";                                            % missing name is 2nd-most common
journals(find(journals == "")) = "Other";                                  % replace missing journals

figure(3);
hold on;
plot([2012,2012],[0,1],'-','Color',[0.6,0.6,0.6],'LineWidth',1);
n=0;
for j = 1:5
    k = mod(j-1,5)+1;
    l = ceil(j./5);
    ix = find(strcmp(journals,journals_unique{jix(j),1})==1);              % list of articles from this journal
    plot(d(ix,3),d(ix,1),[plotcols{k},'o'],'MarkerFaceColor',plotcols{k},'MarkerSize',5);
    mdl = fitlm(d(ix,3),d(ix,1));
    xpred = linspace(min(d(ix,3)),max(d(ix,3)))';
    [ypred,ci] = predict(mdl,xpred);
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


%% REPEAT FOR MEANS OF TOP 25 JOURNALS______________________________________
figure(4);
hold on;
n=0;
for j = 1:25
    k = mod(j-1,5)+1;
    l = ceil(j./5);
    ix = find(strcmp(journals,journals_unique{jix(j),1})==1);               % list of articles from this journal
    M = mean(d(ix,1));
    S = std(d(ix,1));
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


%% COMPARE CLINICAL NEUROPHYSIOLOGY WITH OTHERS_____________________________
% top k journals, with 10 or more papers (1 = Clin Neurophysiol, which published the Chipchase questionnaire)
k = 10;
M = cell2mat(journals_unique(jix(1:k,1),3));                                % means
S = cell2mat(journals_unique(jix(1:k,1),4));                                % SDs
N = cell2mat(journals_unique(jix(1:k,1),2));                                % Ns
V = (N-1).*S.^2;

% comparisons with most popular Journal (#1)
x = 1;                                                                      % comparison = Clin Neurophysiol
y = [2:k];
diff = M(x) - M(y);                                                         % difference from x
Sp = sqrt((V(x) + V(y)) ./ (N(x) + N(y) - 2));                              % pooled SD with x
T = diff ./ (Sp .* sqrt((1./N(x)) + (1./N(y))));                            % unpaired T-test with x
p = 2.*(1-tcdf(abs(T),N(y)-1));                                             % two-tailed p-value with x

% comparisons with best-reportin Journal (#3)
x = 3;                                                                      % comparison = Brain stimulation
y = [1:2,4:k];
diff = M(x) - M(y);                                                         % difference from x
Sp = sqrt((V(x) + V(y)) ./ (N(x) + N(y) - 2));                              % pooled SD with x
T = diff ./ (Sp .* sqrt((1./N(x)) + (1./N(y))));                            % unpaired T-test with x
p = 2.* (1-(tcdf(abs(T),N(1) + N(y)-2)));                                   % two-tailed p-value with x


%% COMPARE CLINICAL NEUROPHYSIOLOGY WITH PRE-CHIPCHASE REPOTING_____________
CN = find(journals=="Clin Neurophysiol");
pre = d(CN,3)>=2007 & d(CN,3)<2012;                                         % pre = 2007 - 2011, 5 years = 12 papers
post = d(CN,3)>=2013;                                                       % post = 2013 - max = 2017, 5 years = 15 papers
M1 = mean(d(CN(pre),1));
M2 = mean(d(CN(post),1));
S1 = std(d(CN(pre),1));
S2 = std(d(CN(post),1));
N1 = sum(pre);
N2 = sum(post);

% unpaired t-test pre vs post publication of the Chipchase (2012) checklist
diff = M2 - M1;
Sp = sqrt((((N1-1).*S1.^2) + ((N2-1).*S2.^2)) ./ (N1 + N2 - 2));
T = diff ./ (Sp.*sqrt( (1./N1) + (1./N2)));
p = 2.* (1-(tcdf(abs(T),N1+N2-2)));