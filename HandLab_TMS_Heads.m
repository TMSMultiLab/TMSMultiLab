%% read the TMSHeads table dumped from LabMan mysql database (participant database used by The Hand Lab)
heads=readtable('HandLab_TMSHeads.csv');

% remove non-head measurements
idx=~strcmp(heads.headtype,'head');
heads(idx,:)=[];

% find unique participant numbers
ps=unique(heads.participantid);

% matrices for data
headstats=nan(numel(ps),17); % 1=participant id;
                             % 2=N-I mean; 3=N-I SD; 4=N-I min; 5=N-I max; 6=N-I n
                             % 7=E-E mean; 8=E-E SD; 9=E-E min; 10=E-E max; 11=E-E n
                             % 12=N-E-I mean; 13=N-E-I SD; 14=N-E-I min; 15=N-E-I max; 16=N-E-I n
                             % 17=ARM mean; 18=ARM SD; 19=ARM min; 20=ARM max; 21=ARM n

%% extract data per participant
% convert table to array
headdata=str2double(table2array(heads(:,5:9)));
n=0;
for p=1:numel(ps)
    disp([' Participant ',int2str(ps(p))]);
    n=n+1;
    idx=heads.participantid==ps(p);
    headstats(n,1)=ps(p);
    for stat=1:4
        start=(stat-1).*5+1;
        headstats(n,start+1)=nanmean(headdata(idx,stat));
        headstats(n,start+2)=nanstd(headdata(idx,stat));
        headstats(n,start+3)=nanmin(headdata(idx,stat));
        headstats(n,start+4)=nanmax(headdata(idx,stat));
        headstats(n,start+5)=sum(isfinite(headdata(idx,stat)));
    end
end

% histograms of head stats
figure(1);
subplot(3,1,1);
title('Distributions of head measurements');
histogram(headstats(:,2));
axis([22,44,0,80]);
xlabel('Nasion - Inion, cm');
subplot(3,1,2);
histogram(headstats(:,7));
axis([22,44,0,80]);
xlabel('Between Pre-auricular points, cm');
subplot(3,1,3);
histogram(headstats(:,12));
axis([22,44,0,80]);
xlabel('Nasion - Pre-auricular - Inion, cm');
set(gcf,'Position',[0,0,400,800]);
print('TMSMultiLab_Head_measurement_distribution.png','-dpng');

% correlations between measures
figure(1);
hold on;
plot(headstats(:,2),headstats(:,7),'ko');
mdl=fitlm(headstats(:,2),headstats(:,7));
[ypred,ci]=predict(mdl,linspace(nanmin(headstats(:,2)),nanmax(headstats(:,2)))');
plot(linspace(nanmin(headstats(:,2)),nanmax(headstats(:,2)))',ypred,'r-');
plot(linspace(nanmin(headstats(:,2)),nanmax(headstats(:,2)))',ci(:,1),'r--');
plot(linspace(nanmin(headstats(:,2)),nanmax(headstats(:,2)))',ci(:,2),'r--');
xlabel('Nasion - Inion, cm');
ylabel('Between Pre-auricular points, cm');
axis([30,44,30,44
set(gcf,'Position',[0,0,700,700]);
print('TMSMultiLab_Head_measurement_correlation.png','-dpng');

idx_ni=isfinite(headstats(:,2));
idx_ee=isfinite(headstats(:,7));
idx_arm=isfinite(headstats(:,17));
idx=logical(idx_ni.*idx_ee);
corrcoef(headstats(idx,2),headstats(idx,7))
idx=logical(idx_ni.*idx_arm);
corrcoef(headstats(idx,2),headstats(idx,17))
idx=logical(idx_ee.*idx_arm);
corrcoef(headstats(idx,7),headstats(idx,17))