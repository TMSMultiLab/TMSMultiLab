%% Implement the Boyles et al. (2026) latency detection algorithm (https://doi.org/10.1038/s41598-026-61560-0)
% Version 1.0, 14th August 2026, written by Nick Holmes @TheHandLab
% 
% usage:
% [onset, stats, options] = MEP_onset_Boyles_2026(data, samplehz, tmstime, grandmean [, options])
%
% data = required. must be nx1 vector of data
% samplehz = required. must be a number, ideally 2000 or greater, for samples collected per second
% tmstime = required. must be a number, how many milliseconds into the data was the TMS pulse presented?
% grandmean = required. must be the same size as data, and comprise the grand mean across all relevant epochs (from the same muscle, participant, condition, etc)
% options = optional. structure with optional parameters; the defaults according to the original article are:
%    options.searchwindow = [5, 45];                                   % default search window = 5 to 45 ms after the TMS pulse
%    options.blocklength = 5;                                          % how long is each block of data to average the derivitives, in samples
%    options.baseline = [100, 1];                                      % baseline window start and stop relative to TMS, in ms before TMS
%    options.amplitudegate = 1.1;                                      % ratio of baseline:data peak-to-peak differences for MEP detection
%    options.peakjitter = 15;                                          % maximum difference between MEP peak and grand average peak
%    options.peakwindowlength = 1.75;                                  % multiply the peak-to-trough duration by this to set search window
%    options.error = 10.^-6;                                           % minimum possible ratio, to avoid dividing by zero
%    options.ratiocutoff = 0.85;                                       % minimum post:pre derivative ratio relative to the peak ratio
%    options.maxlatency = 35;                                          % maximum possible latency
%    options.derivcheck = 4;                                           % samples after the point to check for large derivatives
%    options.basederivSDs = 1.5;                                       % how many SDs for the baseline derivative cutoff?
%    options.derivcheckwindowlength = 2;                               % multiply the peak-to-trough duration by this to set the derivative check window
%    options.plot = false;                                             % don't plot data by default

function [onset, stats, options] = MEP_onset_Boyles_2026(data, samplehz, tmstime, grandmean, options)

    %% CHECK THE INPUT ARGUMENTS__________________________________________
    
    % check number of arguments passed____________________________________
    if nargin<4
        error('MEP_onset_Boyles_2026: at least 4 arguments required (data, samplehz, tmstime, grandmean');
    end
    
    % check size of data__________________________________________________
    if ndims(data)>2 | min(size(data))>1
        error('MEP_onset_Boyles_2026: data must be an nx1 array of samples');
    elseif size(data,2) > size(data,1)
        data = data';
    end
    
    % check sampling frequency____________________________________________
    if numel(samplehz)>1 | ~isnumeric(samplehz) | ~isfinite(samplehz)
        error('MEP_onset_Boyles_2026: samplehz must be a single number');
    end

    % check TMS time______________________________________________________
    if numel(tmstime)>1 | ~isnumeric(tmstime) | ~isfinite(tmstime)
        error('MEP_onset_Boyles_2026: tmstime must be a single number');
    elseif tmstime.*(samplehz./1000) > numel(data)
        error('MEP_onset_Boyles_2026: tmstime is after the last sample in the data');
    end

    % check grandmean_____________________________________________________
    if ndims(grandmean)>2 | min(size(grandmean))>1
        error('MEP_onset_Boyles_2026: grandmean must be an nx1 array of samples');
    elseif size(grandmean,2) > size(grandmean,1)
        grandmean = grandmean';
    end
    if numel(data) ~= numel(grandmean)
        error('MEP_onset_Boyles_2026: grandmean must be the same length as the data');
    end


    %% PROCESS THE OPTIONS________________________________________________
    if nargin==4                                                          % if no options are provided, use defaults:
        options.searchwindow = [5, 45];                                   % default search window = 5 to 45 ms
        options.blocklength = 5;                                          % how long is each block of data to average the derivitives, in samples
        options.baseline = [100, 1];                                      % baseline window start and stop relative to TMS, in ms before TMS
        options.amplitudegate = 1.1;                                      % ratio of baseline:data peak-to-peak differences for MEP detection
        options.peakjitter = 15;                                          % maximum difference between MEP peak and grand average peak
        options.peakwindowlength = 1.75;                                  % multiply the peak-to-trough duration by this to set search window
        options.error = 10.^-6;                                           % minimum possible ratio, to avoid dividing by zero
        options.ratiocutoff = 0.85;                                       % minimum post:pre derivative ratio relative to the peak ratio
        options.maxlatency = 35;                                          % maximum possible latency
        options.derivcheck = 4;                                           % samples after the point to check for large derivatives
        options.basederivSDs = 1.5;                                       % how many SDs for the baseline derivative cutoff?
        options.derivcheckwindowlength = 2;                               % multiply the peak-to-trough duration by this to set the derivative check window
        options.plot = false;                                             % don't plot data by default
    else
        if ~isfield(options, 'searchwindow')
            options.searchwindow = [5, 45];                               % default search window = 5 to 45 ms
        end
        if ~isfield(options, 'blocklength')
            options.blocklength = 5;                                      % how long is each block of data to average the derivitives, in samples
        end
        if ~isfield(options, 'baseline')
            options.baseline = [100, 1];                                  % baseline window start and stop relative to TMS, in ms before TMS
        end
        if ~isfield(options, 'amplitudegate')
            options.amplitudegate = 1.1;                                  % ratio of baseline:data peak-to-peak differences for MEP detection
        end
        if ~isfield(options, 'peakjitter')
            options.peakjitter = 15;                                      % maximum difference between MEP peak and grand average peak
        end
        if ~isfield(options, 'peakwindowlength')
            options.peakwindowlength = 1.75;                              % multiply the peak-to-trough duration by this to set search window
        end
        if ~isfield(options, 'error')
            options.error = 10.^-6;                                       % minimum possible ratio, to avoid dividing by zero
        end
        if ~isfield(options, 'ratiocutoff')
            options.ratiocutoff = 0.85;                                   % minimum post:pre derivative ratio relative to the peak ratio
        end
        if ~isfield(options, 'maxlatency')
            options.maxlatency = 35;                                      % minimum post:pre derivative ratio relative to the peak ratio
        end
        if ~isfield(options, 'derivcheck')
            options.derivcheck = 4;                                       % samples after the point to check for large derivatives
        end
        if ~isfield(options, 'basederivSDs')
            options.basederivSDs = 1.5;                                   % how many SDs for the baseline derivative cutoff?
        end
        if ~isfield(options, 'derivcheckwindowlength')
            options.derivcheckwindowlength = 2;                           % multiply the peak-to-trough duration by this to set the derivative check window
        end
        if ~isfield(options, 'plot')
            options.plot = false;                                         % don't plot data by default
        end
    end
    

    %% PROCESS THE INPUT ARGUMENTS________________________________________
    samples = size(data,1);                                               % samples of data
    tmstimes = round(tmstime .* (samplehz ./ 1000));                      % convert tmstime to discrete samples
    blocklength = round(5.*(samplehz./1000));                             % block length in samples - for slope window
    start = tmstimes + round(options.searchwindow(1) .* samplehz./1000);  % start of search window, in samples
    finish = tmstimes + round(options.searchwindow(2) .* samplehz./1000); % end of search window, in samples
    base.on = (tmstime-options.baseline(1)) .* (samplehz./1000);          % baseline EMG window start, in samples
    base.off = (tmstime-options.baseline(2)) .* (samplehz./1000);         % baseline EMG window stop, in samples


    %% PROCESS THE DATA___________________________________________________
        
    % 0: SUBTRACT BASELINES_______________________________________________
    data = data - nanmean(data(base.on:base.off));                        % subtract mean baseline from data
    grandmean = grandmean - nanmean(grandmean(base.on:base.off));         % subtract mean baseline from grand mean
    
    
    %% 1: DETECT PEAKS____________________________________________________
    
    % in the data_________________________________________________________
    [amp1, lat1] = max(data(start:finish));                               % max positive peak in data
    [amp2, lat2] = min(data(start:finish));                               % max negative peak in data
    [peakamp, peaklat] = max(abs(data(start:finish)));                    % max overall peak in data
    peaklat = peaklat + start;                                            % correct for start of window
    p2p = amp1 - amp2;                                                    % the maximum peak-to-peak in the data

    % in the grand mean___________________________________________________
    [amp3, lat3] = max(grandmean(start:finish));                          % maximum potential MEP peak in grand mean
    [amp4, lat4] = min(grandmean(start:finish));                          % minimum potential MEP peak in grand mean
    grandp2p = amp3-amp4;                                                 % potential MEP peak-to-peak in grand mean
    [grandamp, grandlat] = max(abs(grandmean(start:finish)));             % max overall peak in grandmean data
	grandlat = grantlat + start;                                          % correct for start of window
    
    % in the baseline_____________________________________________________
    base.p2p = max(data(base.on:base.off)) - min(data(base.on:base.off));  % the peak-to-peak in the baseline
    base.derivmean = mean(abs(diff(data(base.on:base.off))));              % baseline mean derivative
    base.derivSD = std(abs(diff(data(base.on:base.off))));                 % baseline SD derivative
    base.derivcutoff = base.derivmean + options.basederivSDs.*base.derivSD;% cutoff for baseline derivatives

    
    %% 2: SET THE AMPLITUDE GATE__________________________________________
    % p2p amplitude in search window post-TMS is > amptlitudegate X the p2p of baseline
    % peak occurs within 15ms of the grand-average peak across all epochs per condition and muscle (and participant?)
    ampgate = options.amplitudegate .* base.p2p;                           % cutoff for amplitude gate
    
    % apply  amplitude gate                     & peak jitter criterion___
    if p2p > ampgate & abs(peaklat - grandlat) <= (options.peakjitter .* (samplehz./1000))
                
        %% MEP DETECTED, SO CONTINUE PROCESSING___________________________
        
        %% 3: DETERMINE SEARCH WINDOWS____________________________________

        % anchor = first peak_____________________________________________
        anchor = min(lat1, lat2) + start;                                 % samples to the first peak in the data
        
        % deltap2p = time between max and min peaks in data_______________
        deltap2p = abs(lat1 - lat2);                                      % samples between detected peaks
        
        % kmin = earliest point to start searching from___________________
        kmin = max(anchor - options.peakwindowlength.*deltap2p, start);   % the latest sample of the peak - search window OR the earliest-possible start
                    
        
        %% 4: FIND THE EARLIEST POINT WHERE THE SIGNAL INCREASES STEEPLY__
        R = nan(samples,3);                                               % variable for storing derivative ratios; N derivatives that are above baseline; mean derivative after sample
        
        % get ratio of slopes after:before each posssible point___________
        for l = anchor : -1 : kmin                                        % go back from the first peak to the minimum sensible latency
            mufwd = mean(abs(diff(data(l+1:l+options.blocklength))));     % mean derivitive in front of sample
            mubwd = mean(abs(diff(data(l-1:-1:l-options.blocklength))));  % mean derivitive behind the sample
            R(l,1) = mufwd ./ (mubwd+options.error);                      % ratio of mean derivatives
            R(l,2) = sum(abs(diff(data(l+1:l+options.derivcheck)))>base.derivmean);% how many of these derivatives exceed the baseline mean?
            R(l,3) = mean(abs(diff(data(l+1:l+options.derivcheckwindowlength.*deltap2p))));% mean derivatives after the sample
        end
                    
        % if slope after timepoint > slope before -> possible onset_______
        % Candidate block
        [Rmax, Rix] = nanmax(R(:,1));                                     % find peak ratio
        Rposs = diff(R(:,1)>(options.ratiocutoff.*Rmax));                 % possible samples - those with at least the minimum proportion of the peak ratio
        Rprimary = find(Rposs(1:Rix)==1, 1, 'Last') + 1;                  % first sample that is at least the criterion ratio
        Rprimaryms = (Rprimary ./ (samplehz./1000)) - tmstime;            % convert this to ms relative to TMS


        %% 5: APPLY ADDITIONAL PARAMETERS_________________________________

        % maximum latency_________________________________________________
        if Rprimaryms <= options.maxlatency
            onset = Rprimaryms;                                            % valid onset so far...
        else
            onset = NaN;                                                   % not valid - too late
        end
        
        % Initial slope - Derivatives of at least 3 of the first 4 samples after k exceed the baseline mean derivative for the current epoch (ensures that this captures genuine, consistent slope)
        if R(Rprimary,2)<3
            onset = NaN;                                                   % not valid - derivatives are too low
        end
	
        % Overall slope - For the window after k (width Dptp x 2), check that mean derivative is greater than baseline mean derivative + 1.5 SD (ensures that overall slope of section after onset is much greater than baseline, indicating a likely MEP)
        if R(Rprimary,3)<base.derivcutoff
            onset = NaN;                                                   % not valid - derivatives are too low
        end


        %% PROVIDE OUTPUT VARIABLES_______________________________________
        if nargout>1
            stats.R = R;
            stats.Rcandidate = Rprimary;
            stats.Rcandidatems = Rprimaryms;
            stats.baseline = base;
        end


        %% PLOT THE DATA?_________________________________________________
        if options.plot
            figure;
            subplot(2,1,1);                                               % top = the data
            hold on;
            xrange = -tmstime : (1000./samplehz) : (samples-1)./(samplehz./1000) - tmstime;% the x-axis time labels
            plot([xrange(1),xrange(end)], [0,0], 'k-');                   % baseline mean
            plot(xrange, data, 'b-');                                     % plot the data
            a = axis;
            plot([-options.baseline(1),-options.baseline(1)], [a(3),a(4)], 'k--');% baseline start
            plot([-options.baseline(2),-options.baseline(2)], [a(3),a(4)], 'k--');% baseline end
            plot([options.searchwindow(1),options.searchwindow(1)], [a(3),a(4)], 'r--');% search window start
            plot([options.searchwindow(2),options.searchwindow(2)], [a(3),a(4)], 'r--');% search window end
            plot([onset,onset], [a(3),a(4)], 'b-');                       % onset detected
            ylabel('EMG amplitude, mV');
            axis([-options.baseline(1).*1.05,100,a(3),a(4)]);             % reset axis
	    
            subplot(2,1,2);                                               % bottom = the diagnostics
            hold on;
            plot(xrange, R(:,1), 'k-');                                   % the derivative ratios tested
            plot([xrange(1),xrange(end)], [base.derivmean,base.derivmean], 'k--');% baseline mean derivative
            set(gca, 'YScale', 'log');
            a = axis;
            axis([-options.baseline(1).*1.05,100,a(3),a(4)]);             % reset axis
            xlabel('Time after TMS, ms');
            ylabel('Derivative ratio, A.U.');
        end
end
