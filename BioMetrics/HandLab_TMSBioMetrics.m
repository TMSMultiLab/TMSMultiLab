%% read the TMSHeads table dumped from LabMan mysql database (participant database used by The Hand Lab)

%% SQL statement: 
% "SELECT headid,headtype,LabMan.TMSHeads.participantid,headdate,IFNULL(nasioninion,\"\") AS \"nasioninion\",IFNULL(intertragal,\"\") AS \"intertragal\",IFNULL(nasionearinion,\"\") AS \"nasionearinion\",IFNULL(LabMan.TMSHeads.armlength,\"\") AS \"armlength\",IFNULL(wristcirc,\"\") AS \"wristcirc\",IFNULL(LabMan.Participants.armlength,\"\") AS \"P_armlength\",IFNULL(height,\"\") AS \"height\",IFNULL(weight,\"\") AS \"weight\",IFNULL(ethnicity,\"\") AS \"ethnicity\"
% FROM LabMan.TMSHeads LEFT JOIN LabMan.Participants ON LabMan.TMSHeads.participantid=LabMan.Participants.participantid
% INTO OUTFILE '/var/www/html/upload/mysql/HandLab_TMSHeads.csv' FIELDS TERMINATED BY ',' ENCLOSED BY '\"' LINES TERMINATED BY '\n'"

%% load the data
heads=readtable('data/HandLab_TMSHeads.csv');

%% rename the variables
heads=renamevars(heads,["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14"],["headid","headtype","participantid","sex","ethnicity","age","nasioninion","intertragal","nasionearinion","armlength","wristcirc","P_armlength","height","weight"]);

% remove non-head measurements
idx=~strcmp(heads.headtype,'head');
heads(idx,:)=[];

% find unique participant numbers
ps=unique(heads.participantid);
ps=ps(isfinite(ps));

% matrices for data
headstats=nan(numel(ps),17); % 1=participant id;
                             %1 6:   2=age mean;    3=age SD;    4=age min;    5=age max;    6=age n
                             %2 7:   7=N-I mean;    8=N-I SD;    9=N-I min;   10=N-I max;   11=N-I n
                             %3 8:  12=E-E mean;   13=E-E SD;   14=E-E min;   15=E-E max;   16=E-E n
                             %4 9:  17=N-E-I mean; 18=N-E-I SD; 19=N-E-I min; 20=N-E-I max; 21=N-E-I n
                             %5 10: 22=ARM mean;   23=ARM SD;   24=ARM min;   25=ARM max;   26=ARM n
                             %6 11: 27=wrist mean; 28=wrist SD; 29=wrist min; 30=wrist max; 31=wrist n
                             %7 12: 32=P_ARM mean; 33=P_ARM SD; 34=P_ARM min; 35=P_ARM max; 36=P_ARM n
                             %8 13: 37=height mean;38=height SD;39=height min;40=height max;41height n
                             %9 14: 42=weight mean;43=weight SD;44=weight min;45=weight max;46=weight n

%% extract data per participant
headdata=table2array(heads(:,6:14));                                                  % convert table to array
n=0;
for p=1:numel(ps)
    disp([' Participant ',int2str(ps(p))]);
    n=n+1;
    idx=heads.participantid==ps(p);
    headstats(n,1)=ps(p);
    for stat=1:9
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
vol=(headstats(:,ni)./pi).*(headstats(:,ee)./pi).*(nanmean(headstats(:,[ni,ee]),2)./pi).*pi.*(2./3000);% estimate half-volume of head, in litres = a*b*c*pi*4/3, where a,b,c, are the three radii
idx_vol=isfinite(vol);
arm=22;
wri=27;
par=32;
idx_par=isfinite(headstats(:,par));                                                   % arm length (participant value)
hei=37;
idx_hei=isfinite(headstats(:,hei));                                                   % height
wei=42;
idx_wei=isfinite(headstats(:,wei));                                                   % weight
jitter=(rand(size(headstats,1),1)-0.5)./3;                                            % jitter up to 0.16cm either way

% histograms of BioMetrics
figure(1);
subplot(2,2,1);
title('Distributions of head measurements');
histogram(headstats(:,ni),20);
axis([28,44,0,80]);
xlabel('Nasion - Inion, cm');

subplot(2,2,2);
histogram(headstats(:,ee),20);
axis([28,44,0,80]);
xlabel('Between Pre-auricular points, cm');

subplot(2,2,3);
histogram(headstats(:,hei),20);
axis([140,200,0,10]);
xlabel('Height, cm');

subplot(2,2,4);
histogram(headstats(:,par),20);
axis([50,100,0,25]);
xlabel('Arm length, cm');

set(gcf,'Position',[0,0,800,800]);
print('data/HandLab_TMSBioMetrics_Distributions.png','-dpng');
close(1);

% correlations between main measures
figure(1);
subplot(2,2,1);                                                                       % NI & EE
hold on;
plot(headstats(:,ni)+jitter(randperm(size(jitter,1))),headstats(:,ee)+jitter(randperm(size(jitter,1))),'k+');
mdl=fitlm(headstats(:,ni),headstats(:,ee));
[ypred,ci]=predict(mdl,linspace(nanmin(headstats(:,ni)),nanmax(headstats(:,ni)))');
plot(linspace(nanmin(headstats(:,ni)),nanmax(headstats(:,ni)))',ypred,'r-');
plot(linspace(nanmin(headstats(:,ni)),nanmax(headstats(:,ni)))',ci(:,1),'r--');
plot(linspace(nanmin(headstats(:,ni)),nanmax(headstats(:,ni)))',ci(:,2),'r--');
xlabel('Nasion - Inion, cm');
ylabel('Between Pre-auricular points, cm');
idx=logical(idx_ni.*idx_ee);                                                          % Ps with both NI and EE measures
r=corrcoef(headstats(idx,ni),headstats(idx,ee));
text(32,43,['r(',int2str(sum(idx)-2),')=',num2str(r(1,2),3)]);
if mdl.Coefficients.Estimate(1)<0
    sign_lbl=' ';
else
    sign_lbl='+';
end
text(32,42,['y=',num2str(mdl.Coefficients.Estimate(2),3),'x ',sign_lbl,num2str(mdl.Coefficients.Estimate(1),3),' cm']);
axis([30,44,30,44]);

subplot(2,2,2);                                                                       % height & weight
hold on;
plot(headstats(:,hei)+jitter(randperm(size(jitter,1))),headstats(:,wei)+jitter(randperm(size(jitter,1))),'k+');
mdl=fitlm(headstats(:,hei),headstats(:,wei));
[ypred,ci]=predict(mdl,linspace(nanmin(headstats(:,hei)),nanmax(headstats(:,hei)))');
plot(linspace(nanmin(headstats(:,hei)),nanmax(headstats(:,hei)))',ypred,'r-');
plot(linspace(nanmin(headstats(:,hei)),nanmax(headstats(:,hei)))',ci(:,1),'r--');
plot(linspace(nanmin(headstats(:,hei)),nanmax(headstats(:,hei)))',ci(:,2),'r--');
xlabel('Height, cm');
ylabel('Weight, kg');
idx=logical(idx_hei.*idx_wei);                                                        % Ps with both measures
r=corrcoef(headstats(idx,hei),headstats(idx,wei));
text(150,105,['r(',int2str(sum(idx)-2),')=',num2str(r(1,2),3)]);
if mdl.Coefficients.Estimate(1)<0
    sign_lbl=' ';
else
    sign_lbl='+';
end
text(150,100,['y=',num2str(mdl.Coefficients.Estimate(2),3),'x ',sign_lbl,num2str(mdl.Coefficients.Estimate(1),3),' kg']);
axis([140,190,40,110]);

subplot(2,2,3);                                                                       % height & arm length
hold on;
plot(headstats(:,hei)+jitter(randperm(size(jitter,1))),headstats(:,par)+jitter(randperm(size(jitter,1))),'k+');
mdl=fitlm(headstats(:,hei),headstats(:,par));
[ypred,ci]=predict(mdl,linspace(nanmin(headstats(:,hei)),nanmax(headstats(:,hei)))');
plot(linspace(nanmin(headstats(:,hei)),nanmax(headstats(:,hei)))',ypred,'r-');
plot(linspace(nanmin(headstats(:,hei)),nanmax(headstats(:,hei)))',ci(:,1),'r--');
plot(linspace(nanmin(headstats(:,hei)),nanmax(headstats(:,hei)))',ci(:,2),'r--');
xlabel('Height, cm');
ylabel('Arm length, cm');
idx=logical(idx_hei.*idx_par);                                                        % Ps with both measures
r=corrcoef(headstats(idx,hei),headstats(idx,par));
text(150,87,['r(',int2str(sum(idx)-2),')=',num2str(r(1,2),3)]);
if mdl.Coefficients.Estimate(1)<0
    sign_lbl=' ';
else
    sign_lbl='+';
end
text(150,85,['y=',num2str(mdl.Coefficients.Estimate(2),3),'x ',sign_lbl,num2str(mdl.Coefficients.Estimate(1),3),' cm']);
axis([140,190,60,90]);

subplot(2,2,4);                                                                       % height & head half-volume
hold on;
plot(headstats(:,hei)+jitter(randperm(size(jitter,1))),vol+jitter(randperm(size(jitter,1))),'k+');
mdl=fitlm(headstats(:,hei),vol);
[ypred,ci]=predict(mdl,linspace(nanmin(headstats(:,hei)),nanmax(headstats(:,hei)))');
plot(linspace(nanmin(headstats(:,hei)),nanmax(headstats(:,hei)))',ypred,'r-');
plot(linspace(nanmin(headstats(:,hei)),nanmax(headstats(:,hei)))',ci(:,1),'r--');
plot(linspace(nanmin(headstats(:,hei)),nanmax(headstats(:,hei)))',ci(:,2),'r--');
xlabel('Height, cm');
ylabel('Head volume, litres');
idx=logical(idx_hei.*idx_vol);                                                        % Ps with both measures
r=corrcoef(headstats(idx,hei),vol(idx));
text(150,4.25,['r(',int2str(sum(idx)-2),')=',num2str(r(1,2),3)]);
if mdl.Coefficients.Estimate(1)<0
    sign_lbl=' ';
else
    sign_lbl='+';
end
text(150,4,['y=',num2str(mdl.Coefficients.Estimate(2),3),'x ',sign_lbl,num2str(mdl.Coefficients.Estimate(1),3),' L']);
axis([140,190,2,4.5]);
set(gcf,'Position',[0,0,800,800]);
print('data/HandLab_TMSBioMetrics_Correlations.png','-dpng');
close(1);

%% RUN THE TMSSites ANALYSIS TOO
HandLab_TMSSites;