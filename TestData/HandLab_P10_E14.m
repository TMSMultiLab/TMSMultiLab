%% Download and process the HandLab dataset P10_E14 for inclusion into the TMSMultiLab TestData set

%% ADD LOCAL TMSMULTILAB REPOSITORY_________________________________________
addpath(genpath('/var/www/html/TMSMultiLab'));

%% RETRIEVE DATA FROM LOCAL DIRECTORY_______________________________________
load('/var/www/html/HandLab/P10_Imitation/P10_E14_BrainConnectivity/data/group/P10_E14_group_data.mat');

% - - - OR - - -

%% DOWNLOAD DATA____________________________________________________________
% https://osf.io/2xytm/files/a2j8e
% Load downloaded data: load('P10_E14_group_data.mat');


%% REMOVE UNNEEDED VARIABLES________________________________________________
clear CI M MEPmax MEPmin N N2 N1 Ns S SE X XCI XS XSE ans autolatency b bs c cb ci;
clear cols colsnum d d1 d2 data diffs ds e f folder fonts fs hl_dir i id ix j;
clear latency latnum lines m markers maxmep mep mpos mx o os p po r s symbols t;
clear trials trials2 trig trigchan triglev tmss v ylimits xrange xrange2;


%% DESCRIPTION OF EXISTING DATASET__________________________________________
% MAIN VARIABLES:
% ERP = 7-dimensional matrix: (Samples x1600, Trials x20, Muscles x4 (FDI-R, ADM-R, FDS-R, EDC-R), Site x2 (M1, IPL), Intensities x7, Conditions x2 (rest, pegboard), Subjects x12);

% MEPs = 8-dimensional matrix: (7 statistics, 2 epochs (post-TMS, pre-TMS), Trials+1, Muscle x4, Site x2, Intensity x7, Condition (rest, pegboard), Subjects x12)
% MEP stats:	% 1 = latency (Bigoni method, but not reliable)
                % 2 = amplitude of the minimal EMG signal (min amp)
		% 3 = latency of the minimal EMG signal (min lat)
		% 4 = amplitude of the maximal EMG signal (max amp)
		% 5 = latency of the maximal EMG signal (max lat)
		% 6 = peak-to-peak MEP amplitude (max-min)
		% 7 = RMS baseline (from -200 to -50ms - before TMS)

% output = 7-dimensional matrix storing the Latencies: (Intensities x7, Epochs (post, pre), Muscles x4, Conditions (rest, pegboard), Sites (M1, SMG), Measures x12, Subjects x12
% latencies were measured on the AVERAGE MEP across 20 trials within each condition -> all EMG data averaged, then MEP measured
% 
% 12 measures:  % 1 = intensity
                % 2 = max mV
                % 3 = max ms
                % 4 = min mV
                % 5 = min ms
                % 6 = latency human 1
                % 7 = latency human 22
                % 8 = auto-onsets: Bigoni et al (2022)
		% 9 = 
                %10 = 
                %11 = 
                %12 = 


%% RE-FORMAT DATA FOR NEW STRUCTURE__________________________________________
% DATA: 2D array:       26880 trials x 1600 samples per trial
% META-DATA: 2D array: 	26880 trials x 16 meta-data variables (1 datasetID, 2 trial, 3 muscle, 4 site, 5 intensity, 6 condition, 7 subject; 8 min-amp, 9 min-lat, 10 max-amp, 11 max-lat, 12 p2p, 13 rms, 14 latency1, 15 latency2, 16 latencyB)

% SUBJECTS: 2D array:   12 x X subject variables (Subject ID, Age, Sex, Handedness, Height, Weight, RMT, AMT,...)

tmp = ERP(:,1:20,:,:,:,:,:);
tmp = reshape(tmp, 1600, 20.*4.*2.*7.*2.*12);
data = tmp';
trials = size(data,1);
meta = nan(trials, 16);


%% FILL IN THE META-DATA IDENTIFIERS________________________
n = 0;
for p = 1:12                                                % for each of 12 participants
    for c = 1:2                                             % and each of 2 conditions
        for I = 1:7	                                    % and each of 7 intensities
            for s = 1:2                                     % and each of two TMS sites
                for m = 1:4                                 % and each of four muscles
                    for t = 1:20                            % and for each of 20 trials...
                        n = n+1;                            % increment the counter (1... 26880)
                        
                        %% DESIGN VARIABLES_________________
                        meta(n,1) = 1;                      % dataset ID
                        meta(n,2) = t;                      % repetition number
                        meta(n,3) = m;                      % muscle
                        meta(n,4) = s;                      % TMS site
                        meta(n,5) = I;                      % TMS intensity
                        meta(n,6) = c;                      % condition (rest, pegboard)
                        meta(n,7) = p;                      % participant number
                        
                        %% EMG / MEP VARIABLES______________
                        meta(n,8) = MEPs(2,1,t,m,s,I,c,p);  % MEP min peak amplitude
                        meta(n,9) = MEPs(3,1,t,m,s,I,c,p);  % MEP min peak latency
                        meta(n,10) = MEPs(4,1,t,m,s,I,c,p); % MEP max peak amplitude
                        meta(n,11) = MEPs(5,1,t,m,s,I,c,p); % MEP max peak latency
                        meta(n,12) = MEPs(6,1,t,m,s,I,c,p); % MEP amplitude (max-min)
                        meta(n,13) = MEPs(7,1,t,m,s,I,c,p); % EMG RMS (baseline)
                        
                        %% LATENCY VARIABLES________________
                        meta(n,14) = output(I,1,m,c,s,6,p); % Human latency estimate 1
                        meta(n,15) = output(I,1,m,c,s,7,p); % Human latency estimate 2
                        meta(n,16) = output(I,1,m,c,s,8,p); % Bigoni et al algorithm latency estimate
			
                    end
                end
            end
        end
    end
end


%% SAVE DATA TO LOCAL DIRECTORY_____________________________
save('HandLab_P10_E14_data.txt', 'data', '-ascii');         % save main data file
save('HandLab_P10_E14_meta.txt', 'meta', '-ascii');         % save meta data file


%% CREATE METADATA TABLE WITH HEADERS_______________________
reps = meta(:,2);
muscles = ms(meta(:,3));
sites = poss(meta(:,4));
intensities = meta(:,5);
conditions = cs(meta(:,6));
participants = meta(:,7);
metadata = table(meta(:,1), reps, muscles', sites', intensities, conditions', participants, meta(:,8), meta(:,9), meta(:,10), meta(:,11), meta(:,12), meta(:,13), meta(:,14), meta(:,15), meta(:,16));
metadata.Properties.VariableNames = {'Dataset', 'Repetition', 'Muscle', 'TMS Site', 'TMS intensity', 'Condition', 'Participant', 'MEP min amp', 'MEP min lat', 'MEP max amp', 'MEP max lat', 'MEP P2P', 'RMS EMG', 'Latency: Human1', 'Latency: Human2', 'Latency: Bigoni'};
writetable(metadata, 'HandLab_P10_E14_meta.csv', 'FileType', 'text');

%% CREATE PARTICIPANTS TABLE WITH HEADERS___________________
sex = {'M','F'};
sex = sex(fem+1);
handed = {'L', 'A', 'R'};
handed = handed(han+1);
participants = table(hs', age', sex', handed', hed(:,1), hed(:,2), rmt', amt', m1(:,1), m1(:,2), smg(:,1), smg(:,2));
participants.Properties.VariableNames = {'ID', 'Age', 'Sex', 'Handedness', 'Nasion-Inion', 'Inter-tragus', 'RMT', 'AMT', 'Site1 lateral', 'Site1 anterior', 'Site2 lateral', 'Site2 anterior'};
writetable(participants, 'HandLab_P10_E14_subjects.csv', 'FileType', 'text');
