%% read the TMSSites table dumped from LabMan mysql database (participant database used by The Hand Lab)

%% SQL statement: 
% SELECT TMSSites.tmssiteid,TMSHeads.headid,headtype,LabMan.Participants.participantid,sex,IFNULL(ethnicity,\"\") AS \"ethnicity\",DATEDIFF(headdate,dateofbirth) AS age,IFNULL(nasioninion,\"\") AS \"nasioninion\",IFNULL(intertragal,\"\") AS \"intertragal\",IFNULL(nasionearinion,\"\") AS \"nasionearinion\",IFNULL(LabMan.TMSHeads.armlength,\"\") AS \"armlength\",IFNULL(wristcirc,\"\") AS \"wristcirc\",IFNULL(LabMan.Participants.armlength,\"\") AS \"P_armlength\",IFNULL(height,\"\") AS \"height\",IFNULL(weight,\"\") AS \"weight\",tmsreference,tmshemisphere,tmssite,tmssitelateral,tmssiteanterior,IFNULL(tmsmuscle,\"NULL\") AS \"tmsmuscle\",tmsmusclestate,tmsmuscleactivation,tmsthresholdmethod,tmsthresholdcriterion,tmsthreshold
% FROM LabMan.TMSSites
% LEFT JOIN LabMan.TMSHeads ON LabMan.TMSSites.headid=LabMan.TMSHeads.headid
% LEFT JOIN LabMan.Participants ON LabMan.TMSHeads.participantid=LabMan.Participants.participantid 
% LEFT JOIN LabMan.TMSThresholds ON LabMan.TMSSites.tmssiteid=LabMan.TMSThresholds.tmssiteid
% INTO OUTFILE '/var/www/html/upload/mysql/HandLab_TMSSites.csv' 
% FIELDS TERMINATED BY ',' ENCLOSED BY '\"' 
% LINES TERMINATED BY '\n';"

%% MUSCLES & NERVES supplied by each nerve (furthest from head = deeper colours)
% Ulnar (blues): ADM, FDI, FCU, FDM, FPB Dark: [36,113,163],[46,134,193],[23,165,137],[19,141,117] 
colours.FDI=[0,0,256];
colours.ADM=[50,50,256];
colours.FDM=[100,100,256];
colours.FCU=[150,150,256];

labels.ADM=[-7,-1]; % text location in the figure
labels.FDI=[-8,0.5];
labels.FDM=[-7,3];
labels.FCU=[-6,3];

% Median (oranges): APB, FCR, FDP, FDS, OP, TE: [186,74,0],[202,111,30],[214,137,16],[212,172,13] 
colours.APB=[256,0,0];
colours.FPB=[256,25,25];
colours.OP=[256,50,50];
colours.TE=[256,75,75];
colours.FDS=[256,100,100];
colours.FCR=[256,150,150];
colours.FDP=[256,200,200];

labels.APB=[-7,1]; % text location in the figure
labels.TE=[-6,2];
labels.FPB=[-6,-0.5];
labels.OP=[-5,3];
labels.FDS=[-6.5,1.5];
labels.FCR=[-3,3];
labels.FDP=[-2,3];

% Radial (purples): BR, ECR, ECU, EDC,TB: [169,50,38],[203,67,53],[136,78,160],[125,60,152]
colours.BR=[256,0,256];
colours.EDC=[256,25,256];
colours.ECR=[256,100,256];
colours.ECU=[256,150,256];
colours.TB=[256,200,256];

labels.BR=[-3,1];
labels.EDC=[-7,0.25];
labels.ECR=[1,3];
labels.ECU=[2,3];
labels.TB=[-4,2];

% Musculocutaneous: BB, Brachialis; greens: [34,153,84],[40,180,99]
colours.BB=[0,256,0];
colours.Brachialis=[25,256,25];

labels.BB=[4,3];
labels.Brachialis=[5,3];

% Axillary & Pectoral: DA,DM (greys) [39,55,70],[46,64,83],[112,123,124],[131,145,146] 
colours.DA=[50,50,50];
colours.DM=[100,100,100];
colours.PM=[150,150,150];   

labels.DA=[-2,-0.5];
labels.DM=[-3.5,-0.5];
labels.PM=[-3,1.5];

% 'hand' (black)
colours.Hand=[0,0,0];
labels.Hand=[-6,-2];

%% read the TMSSites table dumped from LabMan mysql database (participant database used by The Hand Lab)
sites=readtable('data/HandLab_TMSSites.csv');

%% rename the variables
sites=renamevars(sites,["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24","Var25","Var26","Var26"],["tmssiteid","headid","headtype","participantid","sex","ethnicity","age","nasioninion","intertragal","nasionearinion","armlength","wristcirc","P_armlength","P_armspan","height","weight","tmsreference","tmshemisphere","tmssite","tmssitelateral","tmssiteanterior","tmsmuscle","tmsmusclestate","tmsmuscleactivation","tmsthresholdmethod","tmsthresholdcriterion","tmsthreshold"]);
%  1: 5 tmssiteid,headid,headtype*,participantid,sex*
%  6:10 ethnicity*,age,nasioninion,intertragal,nasionearinion
% 11:15 armlength,wristcirc,P_armlength,height,weight
% 16:20 tmsreference*,tmshemisphere*,tmssite*,tmssitelateral,tmssiteanterior
% 21:25 tmsmuscle*,tmsmusclestate*,tmsmuscleactivation,tmsthresholdmethod*,tmsthresholdcriterion*
% 26:30 tmsthreshold
% * = text variables

% remove non-head measurements
idx=~strcmp(sites.headtype,'head');
sites(idx,:)=[];

% find unique participant numbers
ps=unique(sites.participantid);
ps=ps(isfinite(ps));

%% extract data per participant

%% NUMERIC DATA
sitedata=table2array(sites(:,[1,2,4,7:16,20:21,24,27]));                    % convert table to array
%  1: 5 siteid,headid,participantid,age,nasion-inion
%  6:10 inter-tragal,nasion-ear-inion,armlength,wristcirc,P_armlength
% 11:15 height,weight,sitelateral,siteanterior,muscleactivation
% 16:20 threshold

%% TEXT / LABEL DATA
site.labels=table2array(sites(:,[5:6,17:19,22:23,25:26]));                  % extract labels into own table
% 1: 5 sex,ethnicity,tmsreference,tmshemisphere,tmssite
% 6:10 tmsmuscle,tmsmusclestate,thresholdmethod,thresholdmethodcriterion
ms=unique(site.labels(:,6));                                                % list of muscles in dataset (queried in this order)
ms(strcmp(ms,'NULL'))=[];                                                   % remove empty values
site.references=site.labels(:,3);                                           % index of references
site.hemispheres=site.labels(:,4);                                          % index of hemispheres tested
site.label=site.labels(:,5);                                                % index of location tested
site.muscles=site.labels(:,6);                                              % index of muscles tested
site.state=site.labels(:,7);                                                % index of muscle state tested
site.method=site.labels(:,8);                                               % index of threshold method tested

% matrices for data
sitestats=nan(size(sites,1),38);
%  1    participantid
%  2    muscle name
%  3    hemisphere (1=left, 2=right)
% 1  4: 8 age (M,SD,min,max,n) = column 4)
% 2  9:13 nasion-inion (M,SD,min,max,n) = column 5
% 3 14:18 inter-tragal (M,SD,min,max,n) = column 6
% 4 19:23 nasion-inter-tragal-inion (M,SD,min,max,n) = column 7
% 5 24:28 armlength (M,SD,min,max,n) = column 8
% 6 29:33 wristcirc (M,SD,min,max,n) = column 9
% 7 34:38 Parmlength (M,SD,min,max,n) = column 10
% 8 39:43 Parmspan (M,SD,min,max,n) = column 11
% 9 44:48 height (M,SD,min,max,n) = column 12
%10 49:53 weight (M,SD,min,max,n) = column 13
%11 54:58 sitelateral (M,SD,min,max,n) = column 14
%12 59:63 siteanterior (M,SD,min,max,n) = column 15
%13 64:68 muscleactivation (M,SD,min,max,n) = column 16
%14 69:73 threshold (M,SD,min,max,n) = column 17

%% aggregate data per participant__________________________________________
n=0;
for p=1:numel(ps)
    disp([' Participant ',int2str(ps(p)),'...']);
    for m=1:numel(ms)
        %disp(['  Muscle ',int2str(m),': ',ms{m},'...']);
        for h=1:2
            if h==1
                hem='Left';
            else
                hem='Right';
            end
            idx=sitedata(:,3)==ps(p) & strcmp(site.muscles,ms{m}) & strcmp(site.hemispheres,hem) & strcmp(site.references,'Vertex');
            if sum(idx)>0
                n=n+1; 
                sitestats(n,1)=ps(p);                                       % save participant number
                sitestats(n,2)=m;                                           % save muscle number
                sitestats(n,3)=h;                                           % save hemisphere (1=left, 2=right)
                for stat=1:14                                               % for each of the 13 stats in the above table
                    start=3+(stat-1).*5;                                    % starting column in the array (5 measures per stat)
                    sitestats(n,start+1)=nanmean(sitedata(idx,stat+3));     % get the mean
                    sitestats(n,start+2)=nanstd(sitedata(idx,stat+3));      % the SD
                    sitestats(n,start+3)=nanmin(sitedata(idx,stat+3));      % the minimum
                    sitestats(n,start+4)=nanmax(sitedata(idx,stat+3));      % the maximum
                    sitestats(n,start+5)=sum(isfinite(sitedata(idx,stat+3)));% and the number of datapoints averaged
                end
            end
        end
    end
end
sitestats=sitestats(1:n,:);

age=4;
ni=9;
ee=14;
ne=19;
arm=24;
wri=29;
par=34;
spa=39
hei=44;
wei=49;
lat=54;
ant=59;
mus=64;
thr=68;

%% PLOT AVERAGE LOCATIONS FOR M1 HOTSPOTS IN CM & RELATIVE SPACE____________
for f=1:2
    figure(f);
    hold on;
    grid on;
    switch f
        case 1
            plot([-8,8],[0,0],'k-');
            plot([0,0],[-4,4],'k-');
        case 2
            plot([-.25,.25],[0,0],'k-');
            plot([0,0],[-.1,.1],'k-');
    end
    for m=1:numel(ms)
        disp([' Muscle ',int2str(m),': ',ms{m}]);

        % get colour scheme for this muscle
        g=['colour=colours.',ms{m},';'];
        eval(g);
        % get text label location for this muscle
        g=['label=labels.',ms{m},';'];
        eval(g);
        colour=colour./256;
        for h=1:2
            disp(['  hemisphere ',int2str(h)]);
            idx=sitestats(:,2)==m & sitestats(:,3)==h;
            disp(['   ',int2str(sum(idx)),' datapoints found']);
            if sum(idx)>0
                M.lat=nanmean(sitestats(idx,lat));
                M.ant=nanmean(sitestats(idx,ant));
                if f==2
                    M.lat=M.lat./nanmean(sitestats(idx,ee));                % normalise to head size
                    M.ant=M.ant./nanmean(sitestats(idx,ni));                % normalise to head size
                end
                if sum(idx)<6
                    plot(M.lat,M.ant,'s','color',colour,'MarkerSize',12,'LineWidth',1.5);
                end
                if sum(idx)>1
                    SD.lat=nanstd(sitestats(idx,lat));
                    SD.ant=nanmean(sitestats(idx,ant));
                    if f==2
                        SD.lat=SD.lat./nanmean(sitestats(idx,ee));          % normalise to head size
                        SD.ant=SD.ant./nanmean(sitestats(idx,ni));          % normalise to head size
                    end
                    CI.lat=(SD.lat./sum(idx)).*tinv(.975,sum(idx)-1);       % 95% CI
                    CI.ant=(SD.ant./sum(idx)).*tinv(.975,sum(idx)-1);       % 95% CI
                    plot([M.lat-SD.lat,M.lat+SD.lat],[M.ant,M.ant],'-','color',colour,'LineWidth',1.5);
                    plot([M.lat,M.lat],[M.ant-SD.ant,M.ant+SD.ant],'-','color',colour,'LineWidth',1.5);
                    if sum(idx)>5
                        % 95% confidence ellipse
                        t=linspace(0,360,1000);
                        x=CI.lat*sind(t);
                        y=CI.ant*cosd(t);
                        plot(x+M.lat,y+M.ant,'color',colour,'LineWidth',1.5);
                    end
                end
                if h==1
                    if f==1
                        text(label(1),label(2),[ms{m},' (',int2str(sum(idx)),')'],'Color',colour,'FontSize',12);
                    else
                        text(label(1)./36,label(2)./36,[ms{m},' (',int2str(sum(idx)),')'],'Color',colour,'FontSize',12);
                    end
                end
                disp(['  Mean=(',num2str(M.lat,3),',',num2str(M.ant,3),')']);
            end
        end
    end
    set(gcf,'Position',[50,50,1200,600]);
    set(gca,'FontSize',20);
    title('Mean (95% CI) locations of M1 muscles');
    switch f
        case 1
            axis([-9,9,-3.6,3.6]); % mean head size=36cm, this is 25% and 10%
            plot(0,0,'kx','MarkerSize',12);
            text(0.25,-0.25,'Cz','FontSize',12);
            xlabel('Right of vertex, cm');
            ylabel('In front of vertex, cm');
            print('data/HandLab_TMSSites_M1.png','-dpng');
        case 2
            plot(0,0,'kx','MarkerSize',12);
            text(0.005,-0.005,'Cz','FontSize',12);
            plot(-0.2,0,'kx','MarkerSize',12);
            text(-0.205,-0.005,'C3','FontSize',12);
            plot( 0.2,0,'kx','MarkerSize',12);
            text( 0.195,-0.005,'C4','FontSize',12);
            axis([-.25,.25,-.1,.1]);
            xlabel('Right of vertex, relative');
            ylabel('In front of vertex, relative');
            print('data/HandLab_TMSSites_M1_relative.png','-dpng');
    end

end
clear ans colour f h hem idx m n p start stat;

%% histograms of TMS site stats
% figure(1);
% subplot(3,1,1);
% title('Distributions of head measurements');
% histogram(headstats(:,2));
% axis([22,44,0,80]);
% xlabel('Nasion - Inion, cm');
% subplot(3,1,2);
% histogram(headstats(:,7));
% axis([22,44,0,80]);
% xlabel('Between Pre-auricular points, cm');
% subplot(3,1,3);
% histogram(headstats(:,12));
% axis([22,44,0,80]);
% xlabel('Nasion - Pre-auricular - Inion, cm');
% set(gcf,'Position',[0,0,400,800]);
% print('TMSMultiLab_Head_measurement_distribution.png','-dpng');
% 
% % correlations between measures
% figure(1);
% hold on;
% plot(headstats(:,2),headstats(:,7),'ko');
% mdl=fitlm(headstats(:,2),headstats(:,7));
% [ypred,ci]=predict(mdl,linspace(nanmin(headstats(:,2)),nanmax(headstats(:,2)))');
% plot(linspace(nanmin(headstats(:,2)),nanmax(headstats(:,2)))',ypred,'r-');
% plot(linspace(nanmin(headstats(:,2)),nanmax(headstats(:,2)))',ci(:,1),'r--');
% plot(linspace(nanmin(headstats(:,2)),nanmax(headstats(:,2)))',ci(:,2),'r--');
% xlabel('Nasion - Inion, cm');
% ylabel('Between Pre-auricular points, cm');
% axis([30,44,30,44
% set(gcf,'Position',[0,0,700,700]);
% print('TMSMultiLab_Head_measurement_correlation.png','-dpng');
% 
% idx_ni=isfinite(headstats(:,2));
% idx_ee=isfinite(headstats(:,7));
% idx_arm=isfinite(headstats(:,17));
% idx=logical(idx_ni.*idx_ee);
% corrcoef(headstats(idx,2),headstats(idx,7))
% idx=logical(idx_ni.*idx_arm);
% corrcoef(headstats(idx,2),headstats(idx,17))
% idx=logical(idx_ee.*idx_arm);
% corrcoef(headstats(idx,7),headstats(idx,17))