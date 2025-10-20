%% load chipchase reporting data and do analysis on it, to be published here: https://github.com/TMSMultiLab/TMSMultiLab/wiki/Reporting-checklist

% q 1-8     participants
% q 9-25    methods
% q 26-28   paired-pulse
% q 29-30   analysis

%% LOAD THE DATA___________________________________________________________
items_raw=readtable('chipchase_items_raw.csv');
studies=readtable('chipchase_studies.csv');
items_means=readtable('chipchase_items_means.csv');

%% REMOVE DUPLICATE STUDIES________________________________________________
items_raw_duplicates=[];
studies_duplicates=[];

%% ITEMS RAW TABLE_________________________________________________________
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

%% STUDIES TABLE___________________________________________________________
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
jitter=[-0.125,0.125];
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
	s=(q-1).*2+r+4;                                                     % column where the data is
	m=isfinite(table2array(items_means(:,s)));                          % which additional meta-analyses to add?
	n2=nansum(table2array(items_means(m,4)));                           % total N in these meta-analyses
	t2=nansum(round((table2array(items_means(m,4)).*table2array(items_means(m,s)))/100));% total studies meeting the criterion
	
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
text(0.75,1.05,'Participants','FontSize',16);
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

%% TOTAL REPORTING AND CONTROLLING SCORES OVER TIME________________________
% summary data from items_raw table
d=[items_raw.PerR,items_raw.PerC,items_raw.YEAR];                          % all the available data

% add data from studies table
ix=~isfinite(studies.R);                                                   % if 'R' data missing, use TOTAL data
tmp=studies;                                                               % copy temporarily
tmp.R(ix)=tmp.TOTAL(ix);                                                   % copy some TOTAL data into R column
d=[d;tmp.R./100,tmp.C./100,tmp.YEAR];                                      % concatenate all R and T data with items_raw summaries
clear ix tmp;                                                              % remove temporary variables

% index to valid data
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
        yx=logical(ix(:,r).*d(:,3)==y);                                     % index to studies in this year
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