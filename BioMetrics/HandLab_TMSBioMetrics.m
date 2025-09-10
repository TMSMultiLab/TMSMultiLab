%% read the TMSHeads table dumped from LabMan mysql database (participant database used by The Hand Lab)

%% SQL statement: 
% "SELECT headid,headtype,LabMan.TMSHeads.participantid,headdate,IFNULL(nasioninion,\"\") AS \"nasioninion\",IFNULL(intertragal,\"\") AS \"intertragal\",IFNULL(nasionearinion,\"\") AS \"nasionearinion\",IFNULL(LabMan.TMSHeads.armlength,\"\") AS \"armlength\",IFNULL(wristcirc,\"\") AS \"wristcirc\",IFNULL(LabMan.Participants.armlength,\"\") AS \"P_armlength\",IFNULL(LabMan.Participants.armspan,\"\") AS \"P_armspan\",IFNULL(height,\"\") AS \"height\",IFNULL(weight,\"\") AS \"weight\",IFNULL(ethnicity,\"\") AS \"ethnicity\"
% FROM LabMan.TMSHeads LEFT JOIN LabMan.Participants ON LabMan.TMSHeads.participantid=LabMan.Participants.participantid
% INTO OUTFILE '/var/www/html/upload/mysql/HandLab_TMSHeads.csv' FIELDS TERMINATED BY ',' ENCLOSED BY '\"' LINES TERMINATED BY '\n'"

%% load the data
heads=readtable('data/HandLab_TMSHeads.csv');

%% rename the variables
heads=renamevars(heads,["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15"],["headid","headtype","participantid","sex","ethnicity","age","nasioninion","intertragal","nasionearinion","armlength","wristcirc","P_armlength","P_armspan","height","weight"]);

% remove non-head measurements
idx=~strcmp(heads.headtype,'head');
heads(idx,:)=[];

% find unique participant numbers
ps=unique(heads.participantid);
ps=ps(isfinite(ps));

% matrices for data
headstats=nan(numel(ps),17); % 1=participant id;
                             %1   6:  2=age mean;    3=age SD;    4=age min;    5=age max;    6=age n
                             %2   7:  7=N-I mean;    8=N-I SD;    9=N-I min;   10=N-I max;   11=N-I n
                             %3   8: 12=E-E mean;   13=E-E SD;   14=E-E min;   15=E-E max;   16=E-E n
                             %4   9: 17=N-E-I mean; 18=N-E-I SD; 19=N-E-I min; 20=N-E-I max; 21=N-E-I n
                             %5  10: 22=ARM mean;   23=ARM SD;   24=ARM min;   25=ARM max;   26=ARM n
                             %6  11: 27=wrist mean; 28=wrist SD; 29=wrist min; 30=wrist max; 31=wrist n
                             %7  12: 32=P_ARM mean; 33=P_ARM SD; 34=P_ARM min; 35=P_ARM max; 36=P_ARM n
                             %8  13: 37=P_SPAN mean;38=P_SPAN SD;39=P_SPAN min;40=P_SPAN max;41=P_SPAN n
                             %9  14: 42=height mean;43=height SD;44=height min;45=height max;46=height n
                             %10 15: 47=weight mean;48=weight SD;49=weight min;50=weight max;51=weight n

%% extract data per participant
headdata=table2array(heads(:,6:15));                                                  % convert table to array
n=0;
for p=1:numel(ps)
    disp([' Participant ',int2str(ps(p))]);
    n=n+1;
    idx=heads.participantid==ps(p);
    headstats(n,1)=ps(p);
    for stat=1:10
        start=(stat-1).*5+1;
        headstats(n,start+1)=nanmean(headdata(idx,stat));
        headstats(n,start+2)=nanstd(headdata(idx,stat));
        headstats(n,start+3)=nanmin(headdata(idx,stat));
        headstats(n,start+4)=nanmax(headdata(idx,stat));
        headstats(n,start+5)=sum(isfinite(headdata(idx,stat)));
    end
end

%% get indices of useful numbers
age=2;
ni=7;
idx_ni=isfinite(headstats(:,ni));                                                     % Nasion-Inion
ee=12;
idx_ee=isfinite(headstats(:,ee));                                                     % Ear-to-Ear
nei=17;
idx_nei=isfinite(headstats(:,nei));                                                   % Circumference
vol = 52;
headstats(:,vol) = (headstats(:,ni)./pi).*(headstats(:,ee)./pi).*(nanmean(headstats(:,[ni,ee]),2)./pi).*pi.*(2./3000);% estimate half-volume of head, in litres = a*b*c*pi*4/3, where a,b,c, are the three radii
idx_vol=isfinite(headstats(:,vol));
arm=22;
par=32;
idx_par=isfinite(headstats(:,par));                                                   % arm length (participant value)
spa=37;
idx_spa=isfinite(headstats(:,spa));                                                   % arm span (participant value)
hei=42;
idx_hei=isfinite(headstats(:,hei));                                                   % height
wei=47;
idx_wei=isfinite(headstats(:,wei));                                                   % weight
jitter=(rand(size(headstats,1),1)-0.5)./3;                                            % jitter up to 0.16 (cm, L, kg) either way

% histograms of BioMetrics
limits = [5,10,20:20:100,125:25:250,300:50:500];

figure(1);
subplot(2,4,1);
title('Distributions of head measurements');
histogram(headstats(:,ni),20);
a=axis;
N = sum(isfinite(headstats(:,ni)));
b = limits(find(limits>a(4),1)); % first limit-point above current axis (to standardise y-axes)
axis([a(1),a(2),0,b]);
a=axis;
xlabel('Nasion - Inion, cm');
text(a(1) + (a(2)-a(1))/20,b.*.95,['N = ',int2str(N)]);

subplot(2,4,2);
histogram(headstats(:,ee),20);
a=axis;
N = sum(isfinite(headstats(:,ee)));
b = limits(find(limits>a(4),1)); % first limit-point above current axis (to standardise y-axes)
axis([a(1),a(2),0,b]);
a=axis;
xlabel('Between Pre-auricular points, cm');
text(a(1) + (a(2)-a(1))/20,b.*.95,['N = ',int2str(N)]);

subplot(2,4,3);
histogram(headstats(:,nei),20);
a=axis;
N = sum(isfinite(headstats(:,nei)));
b = limits(find(limits>a(4),1)); % first limit-point above current axis (to standardise y-axes)
axis([a(1),a(2),0,b]);
a=axis;
xlabel('Head circumference, cm');
text(a(1) + (a(2)-a(1))/20,b.*.95,['N = ',int2str(N)]);

%subplot(2,4,4);
% histogram(headstats(:,nei),20); - M1-hand location
%a=axis;
%N = sum(isfinite(ni));
%b = limits(find(limits>a(4),1)); % first limit-point above current axis (to standardise y-axes)
%axis([a(1),a(2),0,b]);
%a=axis;
%xlabel('M1-hand lateral, cm');
%text(a(1) + (a(2)-a(1))/20,b.*.95,['N = ',int2str(N)]);

subplot(2,4,5);
histogram(headstats(:,hei),20);
a=axis;
N = sum(isfinite(headstats(:,hei)));
b = limits(find(limits>a(4),1)); % first limit-point above current axis (to standardise y-axes)
axis([a(1),a(2),0,b]);
a=axis;
xlabel('Height, cm');
text(a(1) + (a(2)-a(1))/20,b.*.95,['N = ',int2str(N)]);

subplot(2,4,6);
histogram(headstats(:,wei),20);
a=axis;
N = sum(isfinite(headstats(:,wei)));
b = limits(find(limits>a(4),1)); % first limit-point above current axis (to standardise y-axes)
axis([a(1),a(2),0,b]);
a=axis;
xlabel('Weight, cm');
text(a(1) + (a(2)-a(1))/20,b.*.95,['N = ',int2str(N)]);

subplot(2,4,7);
histogram(headstats(:,par),20);
a=axis;
N = sum(isfinite(headstats(:,par)));
b = limits(find(limits>a(4),1)); % first limit-point above current axis (to standardise y-axes)
axis([a(1),a(2),0,b]);
a=axis;
xlabel('Arm length, cm');
text(a(1) + (a(2)-a(1))/20,b.*.95,['N = ',int2str(N)]);

subplot(2,4,8);
histogram(headstats(:,spa),20);
a=axis;
N = sum(isfinite(headstats(:,spa)));
b = limits(find(limits>a(4),1)); % first limit-point above current axis (to standardise y-axes)
axis([a(1),a(2),0,b]);
a=axis;
xlabel('Arm span, cm');
text(a(1) + (a(2)-a(1))/20,b.*.95,['N = ',int2str(N)]);

set(gcf,'Position',[0,0,1600,800]);
print('data/HandLab_TMSBioMetrics_Distributions.png','-dpng');
close(1);

%% correlations between main measures__________________________________________________
indices = [ ni,  ee;  ni, nei;  ee, nei; hei,vol;  hei, wei; hei,par;  hei,spa;];% pairs for correlation analyses
units =  { 'cm','cm';'cm','cm';'cm','cm';'cm','L';'cm','kg';'cm','cm';'cm','cm';};
subplots = [      1,        2,        3,       4,        5,        6,        7];
lbls = {'Nasion - Inion','Between Pre-auricular points';'Nasion - Inion','Circumference';'Between Pre-auricular points','Circumference';'Height','Head volume';'Height','Weight';'Height','Arm length';'Height', 'Arm span'};
idxs = {logical(idx_ni.*idx_ee);logical(idx_ni.*idx_nei);logical(idx_ee.*idx_nei);logical(idx_hei.*idx_vol);logical(idx_hei.*idx_wei);logical(idx_hei.*idx_par);logical(idx_hei.*idx_spa)};

figure(1);
for s = subplots                                                                       % for each subplot
    subplot(2,4,s);                                                                    % where to plot it
    hold on;
    plot(headstats(:,indices(s,1))+jitter(randperm(size(jitter,1))),headstats(:,indices(s,2))+jitter(randperm(size(jitter,1))),'k+'); % plot the data
    mdl=fitlm(headstats(:,indices(s,1)),headstats(:,indices(s,2)));                    % fit linear model
    [ypred,ci]=predict(mdl,linspace(nanmin(headstats(:,indices(s,1))),nanmax(headstats(:,indices(s,1))))');% get the predicted data and CIs
    plot(linspace(nanmin(headstats(:,indices(s,1))),nanmax(headstats(:,indices(s,1))))',ypred,'r-');% plot the regression line
    plot(linspace(nanmin(headstats(:,indices(s,1))),nanmax(headstats(:,indices(s,1))))',ci(:,1),'r--');% plot the upper CI
    plot(linspace(nanmin(headstats(:,indices(s,1))),nanmax(headstats(:,indices(s,1))))',ci(:,2),'r--');% plot the lower CI
    xlabel([lbls{s,1},', ',units{s,1}]);                                               % x-axis label
    ylabel([lbls{s,2},', ',units{s,2}]);                                               % y-axis label
    r=corrcoef(headstats(idxs{s},indices(s,1)),headstats(idxs{s},indices(s,2)));       % correlation between x and y
    a = axis;                                                                          % get current axis limits
    text(a(1)+(a(2)-a(1))/20,a(4)-(a(4)-a(3))/40,['r(',int2str(sum(idxs{s})-2),')=',num2str(r(1,2),3)]);% add the r-value to the plot
    if mdl.Coefficients.Estimate(1)<0                                                  % get constant sign
        sign_lbl=' ';
    else
        sign_lbl='+';
    end
    text(a(1)+(a(2)-a(1))/20,a(4)-(a(4)-a(3))/10,['y=',num2str(mdl.Coefficients.Estimate(2),3),'x ',sign_lbl,num2str(mdl.Coefficients.Estimate(1),3),' ',units{s,2}]);% add the equation to the plot
    axis(a);                                                                           % set the axis back to where it was
end
set(gcf,'Position',[0,0,1600,800]);                                                    % make figure big
print('data/HandLab_TMSBioMetrics_Correlations.png','-dpng');                          % print figure
close(1);                                                                              % close figure

%% RUN THE TMSSites ANALYSIS TOO
HandLab_TMSSites;
close all;