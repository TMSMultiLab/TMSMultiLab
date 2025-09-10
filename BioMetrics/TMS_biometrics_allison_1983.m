%% Allison et al. (1983) - N=286 study of somatosensory evoked potentials (relevant to TMS-SAI, TMS-PAS)
samplehz=12500;                                                               % data acquired at this frequency

%% PROCESS THE SEP DATA
SEP=readtable('data/allison-t-1983_SEPs.csv');                                % read data from Figures 4, 5, and 6, estimated using G3 data
SEP.ID=NaN(size(SEP,1),1);                                                    % add unique ID column
SEP_labels={'N10','N12a','N12b','N13a','N13b','N14','P15','P16','P18','N20','P22','P26','Height'};

child=SEP.AGE<18;	                                                      % index to children
adult=SEP.AGE>=18;	                                                      % index to adults
male=strcmp(SEP.SEX,'M')==1;                                                  % index to males
female=strcmp(SEP.SEX,'F')==1;                                                % index to females
SEP_idx=nan(size(SEP.Potential,1),numel(SEP_labels)-1);                       % to store all SEP indices

% N10  = ipsilateral shoulder / clavicle, referenced to Fz = median nerve brachial plexus
% N12a-N13a = C7 spinal cord, vs Fz = pre- and post-synaptic spinal afferent volley
% N12b-N13b = C2 spinal cord, vs Fz = pre- and post-synaptic spinal afferent volley
% N14       = C2 spinal cord, vs Fz = more rostral levels of somatosensory afferent pathways
% P15-P18   = P3/P4, vs Fz          = more rostral levels of somatosensory afferent pathways
% N20-P26   = P3/P4, vs Fz          = somatosensory cortex

for p=1:12
    SEP_idx(:,p)=strcmp(SEP.Potential,SEP_labels{p})==1;                              
end

% LOAD REPORTED MEANS
%            Children   Men        Women
%            M    SD    M    SD    M    SD
SEP_means=nan(numel(SEP_labels)-1,3,2);
SEP_means(1,:,:)=[ 8.66,1.10;10.80,0.91; 9.88,0.64]; % N10
SEP_means(2,:,:)=[ 9.85,1.30;12.30,0.86;11.20,0.68]; % N12a
SEP_means(3,:,:)=[10.60,1.30;13.00,0.98;11.80,0.74]; % N12b
SEP_means(4,:,:)=[11.60,1.20;14.40,0.90;13.10,0.71]; % N13a
SEP_means(5,:,:)=[11.90,1.30;14.70,0.95;13.40,0.71]; % N13b
SEP_means(6,:,:)=[12.60,1.20;15.30,0.93;14.00,0.75]; % N14
SEP_means(7,:,:)=[13.40,1.10;16.00,0.99;14.60,0.78]; % P15
SEP_means(8,:,:)=[15.00,1.20;17.50,0.92;16.20,0.79]; % P16
SEP_means(9,:,:)=[16.60,1.30;19.20,0.97;17.90,0.83]; % P18
SEP_means(10,:,:)=[17.50,1.20;20.30,1.10;18.90,0.90];% N20
SEP_means(11,:,:)=[20.70,1.40;23.70,1.50;22.30,1.40];% P22
SEP_means(12,:,:)=[25.10,2.30;27.60,1.80;26.70,1.60];% P26

% ROUND TO NEAREST FEASIBLE NUMBERS
SEP.Latency=round(SEP.Latency.*(samplehz./1000)).*(1000./samplehz);          % round latency data to nearest possible sample
SEP.AGE=round(SEP.AGE);                                                      % histogram of age decimals sugested age was recorded to nearest year, so round

% COMPARE REPORTED AND EXTRACTED MEANS, ADJUST EXTRACTED MEANS TO FIT (most are 0.2-0.3 too high)
for p=1:numel(SEP_labels)-1 % for each labelled potential
    SEP.Latency(child & SEP_idx(:,p))         =SEP.Latency(child & SEP_idx(:,p))         -(mean(SEP.Latency(child & SEP_idx(:,p)))         -SEP_means(p,1,1)); % children
    SEP.Latency(adult & male & SEP_idx(:,p))  =SEP.Latency(adult & male & SEP_idx(:,p))  -(mean(SEP.Latency(adult & male & SEP_idx(:,p)))  -SEP_means(p,2,1)); % adult males
    SEP.Latency(adult & female & SEP_idx(:,p))=SEP.Latency(adult & female & SEP_idx(:,p))-(mean(SEP.Latency(adult & female & SEP_idx(:,p)))-SEP_means(p,2,1)); % adult females
end

%% PROCESS THE HEIGHT DATA
HEIGHT=readtable('data/allison-t-1983_height.csv');                          % read data from Figure 7, estimated using G3 data
HEIGHT.ID=NaN(size(HEIGHT,1),1);                                             % add unique ID column
HEIGHT.Height=HEIGHT.Height-mean(mod(HEIGHT.Height,1));                      % subtract the mean modulus
HEIGHT.Height=round(HEIGHT.Height);                                          % was measured to nearest INCH
HEIGHT.Height=HEIGHT.Height.*2.54;                                           % convert to cm
HEIGHT.AGE=round(HEIGHT.AGE);                                                % histogram of age decimals sugested age was recorded to nearest year, so round
child_height=HEIGHT.AGE<18;	                                             % index to children
adult_height=HEIGHT.AGE>=18;	                                             % index to adults
male_height=strcmp(HEIGHT.SEX,'M')==1;                                       % index to males
female_height=strcmp(HEIGHT.SEX,'F')==1;                                     % index to females
height_means=[146.3,19.8;176.8,7.40;163.1,6.90];                             % reported means in paper, converted to cm
HEIGHT.Height(child_height)                 = HEIGHT.Height(child_height)                 - (mean(HEIGHT.Height(child_height))                 - height_means(1,1)); % children
HEIGHT.Height(adult_height & male_height)   = HEIGHT.Height(adult_height & male_height)   - (mean(HEIGHT.Height(adult_height & male_height))   - height_means(2,1)); % adult males
HEIGHT.Height(adult_height & female_height) = HEIGHT.Height(adult_height & female_height) - (mean(HEIGHT.Height(adult_height & female_height)) - height_means(3,1)); % adult females

%% PROCESS THE DATA
idx=find(ismember(SEP.Potential,'N13a'));                                    % create index for unique participants - match by age and sex; use best-represented data (N10, N13a)
subset=sortrows(SEP(idx,:));                                                 % sorted by sex, age
ID=0;                                                                        % assign unique IDs to unique participants (finds XX unique subjects with YYY datapoints)
for s=1:size(subset,1)                                                       % for each datapoint
    idx=ismember(subset.SEX,subset.SEX(s)) & subset.AGE==subset.AGE(s);      % all datapoints with this sex and age
    if sum(idx)==1                                                           % if only one of this sex and age for this potential
        idx2=ismember(SEP.SEX,subset.SEX(s)) & SEP.AGE==subset.AGE(s);       % index to all datapoints with this sex and age
	if numel(unique(SEP.Potential(idx2))) == numel(SEP.Potential(idx2))  % if only one datapoint for each unique potential
	    ID=ID+1;                                                         % new unique ID
	    SEP.ID(idx2)=ID;                                                 % assign unique ID to all those datapoints

	    % is there a unique height record for this person?
	    idx3=ismember(HEIGHT.SEX,subset.SEX(s)) & HEIGHT.AGE==subset.AGE(s);
	    if numel(unique(HEIGHT.Height(idx3))) == numel(HEIGHT.Height(idx3))
	        HEIGHT.ID(idx3)=ID;
	    end
	end
    end
end

%% plot relationships between age, sex and latency, for each potential
for p=1:numel(SEP_labels)-1
    figure(p);
    for s=1:4
        subplot(2,2,s);
        hold on;
	titletext='';
        if s<3
	    idx=SEP_idx(:,p) & male;
	    plotcol='b';
	    titletext=[titletext,'Male, '];
	else
	    idx=SEP_idx(:,p) & female;
	    plotcol='r';
            xlabel('Age, years');
	    titletext=[titletext,'Female, '];
	end
	if mod(s,2)==1
            ylabel('Latency, ms');
	    X=SEP.AGE(idx);
	    xticks([1,10,20,40,100]);
	    xticklabels({'1','10','20','40','100'});
	    titletext=[titletext,'linear age'];
	else
	    X=log10(SEP.AGE(idx));
	    xticks(log10([0.1,1,5,10,20,40,100]));
	    xticklabels({'0.1','1','5','10','20','40','100'});
	    titletext=[titletext,'log age'];
	end
	Y=SEP.Latency(idx);
        plot(X,Y,[plotcol,'o']);
    
        % fit linear model
        mdl=fitlm(X,Y);
        xpred=linspace(min(X),max(X))';
        [ypred,ci]=predict(mdl,xpred);
        plot(xpred,ypred,[plotcol,'-'],'Linewidth',2);
        plot(xpred,ci(:,1),[plotcol,'--']);
        plot(xpred,ci(:,2),[plotcol,'--']);
	
	a=axis;
	text(mean(a(1:2)),a(3)+(a(4)-a(3))./10,1,['R^2=',num2str(mdl.Rsquared.Adjusted,3)],'Color',plotcol,'FontSize',12);
	title(titletext,'Color',plotcol,'FontSize',10);	
	    
    end   
    print(['data/allison-t_1983-',SEP_labels{p}],'-dpng');
    close(p);
end
    
% 17 unique participants - correlations between all components, N10-P26, with height
idx=isfinite(SEP.ID);
subset=sortrows(SEP(idx,:));
ss=unique(subset.ID);

% build table of potentials for unique subjects
potentials=nan(max(ss),13);
for s=1:max(ss)
    for p=1:12
        idx=find(subset.ID==s & strcmp(subset.Potential,SEP_labels{p}));
        if ~isempty(idx)
            potentials(s,p)=subset.Latency(idx);
        end
    end
    idx=find(HEIGHT.ID==s);
    if ~isempty(idx)
        potentials(s,13)=HEIGHT.Height(idx);
    end
end

% plot it!
xs=[1,2,4,8,10]; % what's on the x-axes?
ys=[13,1,2,4,8]; % what's on the y-axes?
figure(1);
for q=1:5
    for p=1:5
        subplot(5,5,((q-1).*5)+p);
	hold on;
	if xs(p)==ys(q)
	    plotcol='w';
	else
	    plotcol='k';
	end
	X=potentials(:,xs(p));
	Y=potentials(:,ys(q));
        plot(X,Y,[plotcol,'.'],'MarkerSize',8);
	
	% fit linear model
        mdl=fitlm(X,Y);
        xpred=linspace(min(X),max(X))';
        [ypred,ci]=predict(mdl,xpred);
        plot(xpred,ypred,[plotcol,'-']);

	% format graph
        set(gca,'FontSize',6);
	if q==1
            title(SEP_labels{xs(p)});
	elseif q==5
	    xlabel('Latency, ms');
	end
        if p==1
            ylabel(SEP_labels{ys(q)});
        end
	a=axis;
	text(mean(a(1:2)),a(3)+(a(4)-a(3))./10,1,['R^2=',num2str(mdl.Rsquared.Adjusted,3)],'Color',plotcol,'FontSize',6);
    end
end
print('data/allison-t_1983_correls.png','-dpng');
close(1);