%% load & process (TMS) biometrics data

%% Allison et al. (1983) - N=286 study of somatosensory evoked potentials (relevant to TMS-SAI, TMS-PAS)
SEP=readtable('Data/allison-t-1983_SEPS.csv');                                % read data from Figures 4, 5, and 6, estimated using G3 data
height=readtable('Data/allison-t-1983_height.csv');                           % read data from Figure 7, estimated using G3 data
samplehz=12500;                                                               % data acquired at this frequency
data.ID=NaN(size(data,1),1);                                                  % add unique ID column
child=data.AGE<18;	                                                      % index to children
adult=data.AGE>=18;	                                                      % index to adults

% compare mean of extracted data with reported means, and adjust
%        Children   Men        Women
%          M   SD     M   SD     M   SD
heights=[146.3,19.8;176.8,7.40;163.1,6.90];
N10s=   [ 8.66,1.10;10.80,0.91; 9.88,0.64];
N12as=  [ 9.85,1.30;12.30,0.86;11.20,0.68];
N12bs=  [10.60,1.30;13.00,0.98;11.80,0.74];
N13as=  [11.60,1.20;14.40,0.90;13.10,0.71];
N13bs=  [11.90,1.30;14.70,0.95;13.40,0.71];
N14s=   [12.60,1.20;15.30,0.93;14.00,0.75];
P15s=   [13.40,1.10;16.00,0.99;14.60,0.78];
P16s=   [15.00,1.20;17.50,0.92;16.20,0.79];
P18s=   [16.60,1.30;19.20,0.97;17.90,0.83];
N20s=   [17.50,1.20;20.30,1.10;18.90,0.90];
P22s=   [20.70,1.40;23.70,1.50;22.30,1.40];
P26s=   [25.10,2.30;27.60,1.80;26.70,1.60];

data.Latency = round(data.Latency.*(samplehz./1000)).*(1000./samplehz);       % round latency data to nearest possible sample
data.AGE=round(data.AGE);                                                     % histogram of age decimals sugested age was recorded to nearest year, so round
idx=find(ismember(data.Potential,'N13a'));                                    % create index for unique participants - match by age and sex; use best-represented data (N10, N13a)
subset=sortrows(data(idx,:));                                                 % sorted by sex, age
ID=0;                                                                         % assign unique IDs to unique participants (finds XX unique subjects with YYY datapoints)
for s=1:size(subset,1)                                                        % for each datapoint
    idx=ismember(subset.SEX,subset.SEX(s)) & subset.AGE==subset.AGE(s);       % all datapoints with this sex and age
    if sum(idx)==1                                                            % if only one of this sex and age for this potential
        idx2=ismember(data.SEX,subset.SEX(s)) & data.AGE==subset.AGE(s);      % index to all datapoints with this sex and age
	if numel(unique(data.Potential(idx2))) == numel(data.Potential(idx2)) % if only one datapoint for each unique potential
	    ID=ID+1;                                                          % new unique ID
	    data.ID(idx2)=ID;                                                 % assign unique ID to all those datapoints
	end
    end
end