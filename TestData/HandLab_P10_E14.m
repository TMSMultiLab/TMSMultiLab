%% HACKATHON analysis

% download data: % https://osf.io/2xytm/files/a2j8e

% open MATLAB

% Load data load('P10_E14_group_data.mat');

% remove unneeded variables clear CI M MEPmax MEPmin N N2 N1 S SE X XS XSE ans autolatency b bs c cb ci;

clear cols colsnum d d1 d2 data diffs ds e f folder fonts fs hl_dir i id j;

clear latency latnum lines m m1 markers maxmep mep mpos mx o p r s symbols t;

clear trials trials2 trig trigchan triglev tmss v ylimits xrange xrange2;

% MAIN VARIABLES: % ERP = zeros(samples x1600, trials x20, Muscles x4, Site x2 (M1, IPL), Intensities x7, Conditions x2 (rest, active), Subjects x12);% samples, trials, 4 muscles, 2 positions, 7 intensities, 2 conditions(rest/pegboard), Ps)

% MEPs=nan(6, 2, trials+1, muscle, site, Intensity, condition (rest, active), max(ps));% % 7 = 1:latency (Bigoni method, but not reliable) % 2: min amp % 3: min lat % 4: max amp % 5: max lat % 6: MEP amplitude (max-min) % 7: RMS baseline

% OUTPUT = Latencies % 7 intensities, 2 epochs (post, pre), 4 muscles, 2 conditions (rest, pegboard), 2 sites (M1, SMG), 12 measures, 12 participants % 12 measures (intensity; max mV, max ms; min mV, min ms; 6=onset1; 7=onset2; 8=auto-onsets: Bigoni et al (2022)) % 7 2 4 2 2 12 12

%% FORMAT FOR NEW DATA STRUCTURE_______________________________

% NEW DATA... 27000 x 1600; % trials x samples
% META-DATA 27000 x 16 meta-data (1 datasetID, 2 trial, 3 muscle, 4 site, 5 intensity, 6 condition, 7 subject; 8 min-amp, 9 min-lat, 10 max-amp, 11 max-lat, 12 p2p, 13 rms, 14 latency1, 15 latency2, 16 latencyB)

tmp = ERP(:,1:20,:,:,:,:,:);
tmp = reshape(tmp, 1600, 20.*4.*2.*7.*2.*12);
data = tmp';
trials = size(data,1);
meta = nan(trials, 16);

%% fill in the meta-data identifiers
n = 0;
for t = 1:20
    for m = 1:4
        for s = 1:2
            for I = 1:7
	        for c = 1:2
	            for p = 1:12
		        n = n+1;
			
			%% DESIGN VARIABLES_________________
                        meta(n,1) = 1; % dataset ID
			meta(n,2) = t; % repetition number
			meta(n,3) = m; % muscle
			meta(n,4) = s; % TMS site
			meta(n,5) = I; % TMS intensity
			meta(n,6) = c; % condition (rest, pegboard)
			meta(n,7) = p; % participant
			
			%% EMG / MEP variables_____________
			meta(n,8) = MEPs(2,1,t,m,s,I,c,p); % MEP min peak amplitude
			meta(n,9) = MEPs(3,1,t,m,s,I,c,p); % MEP min peak latency
			meta(n,10) = MEPs(4,1,t,m,s,I,c,p);% MEP max peak amplitude
			meta(n,11) = MEPs(5,1,t,m,s,I,c,p);% MEP max peak latency
			meta(n,12) = MEPs(6,1,t,m,s,I,c,p);% MEP amplitude (max-min)
			meta(n,13) = MEPs(7,1,t,m,s,I,c,p);% EMG RMS (baseline)
			
			%% LATENCY VARIABLES_______________
			meta(n,14) = output(I,1,m,c,s,8,p);% Human latency estimate 1
			meta(n,15) = output(I,1,m,c,s,9,p);% Human latency estimate 2
			meta(n,16) = output(I,1,m,c,s,10,p);% Bigoni et al algorithm latency estimate
		    end
		end
	    end
	end
    end
end
