%% EXPORT RAW DATA FROM P10_E14 TO BIDS FORMAT_____________________________
% see README_BIDS.md
% 
addpath(genpath('/var/www/html/HandLabToolBox'));
addpath(genpath('/var/www/html/TMSMultiLab'));
addpath '/var/www/html/HandLab/P10_Imitation/P10_E14_BrainConnectivity';    % for the across-experiment scripts
data_dir = '/var/www/html/TMSMultiLab/TestData';                            % where is the data going to be saved?
hl_dir = 'HandLab_P10_E14';                                                 % top-level study directory
disp('exporting EMG data...');                                              % display info to screen
epoch = 2.0;                                                                % 1000ms before and after TMS
samplehz = 4000;
samples = epoch.*samplehz;
ps =  [       1,       2,      3,        4,        5,      6,        7,      8         9,     10,      11,        12]; % Participant number
hs =  [       1,    1500,   1501,     1502,     1503,   1504,     1505,   1506,     1507,   1508,     1509,     1510]; % subject H numbers
age = [      46,      21,     21,       40,       33,     48,       25,     30,       20,     20,       24,       25]; % age in years
fem = [       0,       0,      0,        1,        0,      1,        1,      1,        0,      1,        0,        0]; % female? 0=male, 1=female
han = [       2,       0,      2,        2,        2,      2,        2,      2,        2,      2,        2,        2]; % left=0, ambi=1, right=2
hed = [ 38,36.5;   36,34;  38,38;    36,34;  40,39.5;  38,36;    36,39;  36,38;    42,40;  38,39;    38,39;    37,36;];% size of head in cm, L-R,N-I
rmt = [      55,      45,     52,       49,       48,     40,       70,    nan,       63,     46,       46,       54]; % resting motor threshold (%MSO)
amt = [      30,      30,     33,       28,       34,     21,       38,     75,       36,     30,       23,       38]; % active motor threshold (%MSO)
m1 =  [-5.5,1.5;-5.5,0.5; -5,0.5;     -5,1;     -5,1; -5.5,0;     -6,1;   -5,2;     -5,1; -5,0.5;     -5,1;     -5,1];% M1 location, (RIGHT,ANTERIOR)
smg=[   -9,-2;   -9,-3;-8.5,-3;-8.5,-2.5;-8.5,-2.5;-9,-3.5;-9.5,-2.5;nan,nan;-8.5,-2.5;-8.5,-3;-8.5,-2.5;-8.5,-2.5];% SMG location, (RIGHT, ANTERIOR)
trials = 20;                                                                % TMS pulses in total, use variable t as index
poss = {'M1','SMG'};                                                        % All TMS positions
fs = [];                                                                    % variable for data files
os = {'1','2','3','4','5','6','7'};                                         % different intensity labels (piggyback on the orientation variable)
ms = {'FDI-R','ADM-R','FDS-R','EDC-R'};                                     % muscles recorded from, use variable m as index
cs = {'Rest','Pegboard'};                                                   % condition label(s)
tmstime = 1000;                                                             % TMS presented at this time, in ms relative to epoch onset
tmstime2 = round(tmstime./1000,2);                                          % TMS time in seconds

% experiment-specific
trigchan = 6;                                                               % number of the trigger channel
triglev = 3;                                                                % trigger level (typically 5V MAX)
emgchans = 2:5;                                                             % where's the EMG in the raw data?
P10_E14_latencies;                                                          % load the (manual) latency data (6=NH; 7=DG)


%% CHECK IF DATA DIRECTORY EXISTS__________________________________________
if exist([data_dir,'/',hl_dir], 'dir')~=7
    mkdir(data_dir,'/',hl_dir);
end


%% CHECK OR CREATE THE FILES_______________________________________________

%%  dataset_description.json_______________________________________________
if exist([data_dir,'/',hl_dir,'/dataset_description.json'], 'file')~=2     % create this manually - easier
    warning('dataset_description.json is missing');
end

%	*_coordsystem.json
%
%	*_electrodes.json
%	*_electrodes.tsv
%
%	*_events.tsv
%	*_events.json
%
%	*_emg.edf
%	*_emg.json


%%  *_channels.tsv - applies to all subjects, so copied into each emg directory
name = {'time', ms{:}, 'TMS'}';
type = {'LATENCY'; 'EMG'; 'EMG'; 'EMG'; 'EMG'; 'TRIG'};
units = {'s'; 'mV'; 'mV'; 'mV'; 'mV'; 'V'};
signal_electrode = {'n/a', ms{:}, 'n/a'}';
reference = {'n/a'; 'bipolar'; 'bipolar'; 'bipolar'; 'bipolar'; 'n/a'};
group = {'n/a'; 'amp1'; 'amp1'; 'amp2'; 'amp2'; 'n/a'};
target_muscle = {'n/a'; 'right FDI'; 'right ADM'; 'right FDS'; 'right EDC'; 'n/a'};
description = {'time to TMS'; 'target muscle'; 'nontarget muscle'; 'nontarget muscle'; 'nontarget muscle'; 'TMS 5V TLL'};
channels = table(name, type, units, signal_electrode, reference, group, target_muscle, description);% build the channels table


%%  participants.tsv & participants.json___________________________________
if exist([data_dir,'/',hl_dir,'/participants.json'], 'file')~=2            % create this manually - easier
    warning('participants.json is missing');
end
participant_id = repmat({''},max(ps),1);                                   % initialise the participant ID variable
species = repmat({'homo sapiens'},max(ps),1);                              % initialise the species variable
sex = {'M','F'};                                                           % codes for 0, 1
sex = sex(fem+1)';                                                         % participant sex
handedness = {'L','A','R'};                                                % codes for 0, 1, 2
handedness = handedness(han+1)';                                           % participant handedness
ear2ear = hed(:,1);                                                        % distance from left to right helix-tragus intersection
nas2in = hed(:,2);                                                         % distance from nasion to inion
m1lat = m1(:,1);
m1ant = m1(:,2);
smglat = smg(:,1);
smgant = smg(:,2);
age = age';
rmt = rmt';
amt = amt';


%% START THE PROCESSING____________________________________________________
for p = ps

    if p==1
        v=2;
    else
        v=1;                                                                % visit =1 for most Ps
    end
    disp(['P',int2str(p),'...']);
    if p<10
        pstr = ['sub-0',int2str(p)];
    else
        pstr = ['sub-',int2str(p)];
    end
    participant_id{p} = pstr;

    % check if the participant directory exists____________________________
    if exist([data_dir,'/',hl_dir,'/',pstr], 'dir')~=7
        mkdir([data_dir,'/',hl_dir,'/',pstr]);
    end
    
    % check if the emg directory exists____________________________________
    if exist([data_dir,'/',hl_dir,'/',pstr,'/emg'], 'dir')~=7
        mkdir([data_dir,'/',hl_dir,'/',pstr,'/emg']);
    end

    % check if the channels file exists____________________________________
    if exist([data_dir,'/',hl_dir,'/',pstr,'/emg/',pstr,'_channels.tsv'], 'file')~=2
        writetable(channels, [data_dir,'/',hl_dir,'/',pstr,'/emg/',pstr,'_channels.tsv'], 'FileType', 'Text', 'Delimiter', '\t');% save channels.tsv
    end
    
    % set up the Scans file________________________________________________
    filename = {''};				                            % filename - with modality folder
    acq_time = {''};				                            % YYYY-MM-DDThh:mm:ss[.000000][Z|+hh:mm|-hh:mm]; YEAR can be shifted to <=1925 to preserve anonymity
    f = 0;                                                                  % for counting filenames

    for c = 1:2                                                             % for each condition (rest, pegboard)
        disp([' C',int2str(c),': ',cs{c},'...']);
        folder = ['data/raw/S',int2str(p),'_H',int2str(hs(p)),'_V',int2str(v),'/'];% folder for the data
        for po = 1:numel(poss)                                              % for each position, M1, SMG
            fs = dir(['/var/www/html/HandLab/P10_Imitation/P10_E14_BrainConnectivity/',folder,'S',int2str(p),'_H',int2str(hs(p)),'*_',poss{po},'_',cs{c},'.txt']);% list the files
            for o = 1:size(fs,1)                                            % for each of the 'orientations' (intensities)
	        f = f+1;                                                    % increment the file counter

                % open the file in .txt format_____________________________
                % data = hl_powerlab_import([fs(o).folder,'/',fs(o).name],7); % import data; 1 file per position, condition, intensity: 1=sample, 2-5=EMG, 6=TMS
		% data = data(:,1:6);                                       % remove redundant trigger channel

                % build scan details_______________________________________
		tmp = fs(o).name;
		filename{f,1} = ['emg/',pstr, '_task-', lower(cs{c}), '_acq-', poss{po}, tmp(7:8), '.edf'];
		acq_time{f,1} = datetime(fs(o).date, 'Format', 'yyyy-MM-dd''T''HH:mm:ss');

                % save the data in BIDS format_____________________________
		
		
                % extract the TMS triggers_________________________________
                %trig = diff(data(:,trigchan)>triglev);                      % samples above trigger level
                %trig(trig==-1) = 0;                                         % remove the ends of triggers
                %tmss = find(diff(trig)==1);                                 % triggers changing from below to above
                %trials2 = size(tmss,1);                                     % trials2 = real number of trials in this dataset

                %for t = 1:trials2                                           % for each of the trials
                    % construct trigger-level event descriptions___________

                %end
            end                                                             % of intensity loop
        end                                                                 % of position loop
    end                                                                     % of condition loop
    
    
    %% WRITE THE SCANS TABLE TO FILE________________________________________
    scans = table(filename, acq_time);                                      % build the scans table
    writetable(scans, [data_dir,'/',hl_dir,'/',pstr,'/',pstr,'_scans.tsv',], 'FileType', 'Text', 'Delimiter', '\t');% save scans.tsv
    
end                                                                         % of subject loop


%% WRITE THE PARTICIPANTS TABLE TO FILE_____________________________________
participants = table(participant_id, species, age, sex, handedness, ear2ear, nas2in, rmt, amt, m1lat, m1ant, smglat, smgant);
writetable(participants, [data_dir,'/',hl_dir,'/participants.tsv'], 'FileType', 'Text', 'Delimiter', '\t');
