%% use the cusum (cumulative sum) method to detect MEP latency
%
% usage:
% [onset, stats, options] = MEP_onset_cusum(data, samplehz, tmstime, [, options])
%
% data = required. must be nx1 vector of data
% samplehz = required. must be a number, ideally 2000 or greater, for samples collected per second
% tmstime = required. must be a number, how many milliseconds into the data was the TMS pulse presented?
% options = optional. structure with optional parameters; the defaults according to the original article are:
%    options.searchwindow = [10, 50];                                  % default search window = 10 to 50 ms after the TMS pulse
%    options.baseline = [100, 1];                                      % baseline window start and stop relative to TMS, in ms before TMS
%    options.climit = 10;                                              % limit in number of SDs from the mean - signal needs to pass this level before detection
%    options.mshift = 3.09;                                            % this number of SDs is subtracted from each sample before adding to the cumulative sum - gate the small fluctuations
%    options.plot = false;                                             % don't plot data by default
%
% Version 1.0, 19th August 2026, written by Nick Holmes @TheHandLab
% cumulative sum is a well-used method in engineering to detect drift or deviation away from the mean in a signal.
% In neurophysiology, it has been used for at least 50 years, for example by: Ellaway (1978):
% Ellaway PH (1978) Cumulative sum technique and its application to the analysis of peristimulus time histograms.
% Electroencephalography and Clinical Neurophysiology, 45(2):302-304 https://doi.org/10.1016/0013-4694(78)90017-2

function [onset, stats, options] = MEP_onset_cusum(data, samplehz, tmstime, options)

    %% CHECK THE INPUT ARGUMENTS__________________________________________
    
    % check number of arguments passed____________________________________
    if nargin<3
        error('MEP_onset_cusum: at least 3 arguments required (data, samplehz, tmstime)');
    end
    
    % check size of data__________________________________________________
    if ndims(data)>2 | min(size(data))>1
        error('MEP_onset_cusum: data must be an nx1 array of samples');
    elseif size(data,2) > size(data,1)
        data = data';
    end
    
    % check sampling frequency____________________________________________
    if numel(samplehz)>1 | ~isnumeric(samplehz) | ~isfinite(samplehz)
        error('MEP_onset_cusum: samplehz must be a single number');
    end

    % check TMS time______________________________________________________
    if numel(tmstime)>1 | ~isnumeric(tmstime) | ~isfinite(tmstime)
        error('MEP_onset_cusum: tmstime must be a single number');
    elseif tmstime.*(samplehz./1000) > numel(data)
        error('MEP_onset_cusum: tmstime is after the last sample in the data');
    end


    %% PROCESS THE OPTIONS________________________________________________
    if nargin==3                                                          % if no options are provided, use defaults:
        options.searchwindow = [10, 50];                                  % default search window = 10 to 50 ms
        options.baseline = [100, 1];                                      % baseline window start and stop relative to TMS, in ms before TMS
        options.climit = 10;                                              % limit in number of SDs from the mean - signal needs to pass this level before detection
        options.mshift = 3.09;                                            % this number of SDs is subtracted from each sample before adding to the cumulative sum - gate the small fluctuations
        options.plot = false;                                             % don't plot data by default
    else
        if ~isfield(options, 'searchwindow')
            options.searchwindow = [10, 50];                              % default search window = 10 to 50 ms
        end
        if ~isfield(options, 'baseline')
            options.baseline = [100, 1];                                  % baseline window start and stop relative to TMS, in ms before TMS
        end
        if ~isfield(options, 'climit')
            options.climit = 10;                                          % limit in number of SDs from the mean - signal needs to pass this level before detection
        end
        if ~isfield(options, 'mshift')
            options.mshift = 3.09;                                        % this number of SDs is subtracted from each sample before adding to the cumulative sum - gate the small fluctuations
        end
        if ~isfield(options, 'plot')
            options.plot = false;                                         % don't plot data by default
        end
    end


    %% PROCESS THE INPUT ARGUMENTS________________________________________
    samples = size(data,1);                                               % samples of data
    tmstimes = round(tmstime .* (samplehz ./ 1000));                      % convert tmstime to discrete samples
    start = tmstimes + round(options.searchwindow(1) .* samplehz./1000);  % start of search window, in samples
    finish = tmstimes + round(options.searchwindow(2) .* samplehz./1000); % end of search window, in samples
    base.on = (tmstime-options.baseline(1)) .* (samplehz./1000);          % baseline EMG window start, in samples
    base.off = (tmstime-options.baseline(2)) .* (samplehz./1000);         % baseline EMG window stop, in samples


    %% PROCESS THE DATA___________________________________________________
        
    % 0: SUBTRACT BASELINES_______________________________________________
    data = data - nanmean(data(base.on:base.off));                        % subtract mean baseline from data
    
    
    % 0: GET BASELINE MEAN AND SD_________________________________________
    base.mean = nanmean(data(base.on:base.off));                          % calculate mean signal over baseline
    base.SD = nanstd(data(base.on:base.off));                             % calculate SD signal over baseline
    
    
    %% 1: RUN CUSUM_______________________________________________________
    [stats.uppers stats.lowers stats.uppersum stats.lowersum] = cusum(data(start:finish), options.climit, options.mshift, base.mean, base.SD);
    % uppper = first index of upper boundary detected, in samples
    % lower  = first index of upper boundary detected, in samples
    % uppersum = upper sum to the boundary
    % lowersum = lower sum to the boundary
    stats.upper = (stats.uppers + start - tmstimes) .* (1000./samplehz);   % convert to ms
    stats.lower = (stats.lowers + start - tmstimes) .* (1000./samplehz);   % convert to ms    
    onset = nanmin([stats.upper,stats.lower]);                             % output the first sample


    %% PLOT THE DATA?_____________________________________________________
    if options.plot
        figure;
        subplot(2,1,1);                                                   % top = the data
        hold on;
        xrange = -tmstime : (1000./samplehz) : (samples-1)./(samplehz./1000) - tmstime;% the x-axis time labels
        plot([xrange(1),xrange(end)], [0,0], 'k-');                       % baseline mean
        plot(xrange, data, 'b-');                                         % plot the data
	a = axis;
	plot([-options.baseline(1),-options.baseline(1)], [a(3),a(4)], 'k--');% baseline start
	plot([-options.baseline(2),-options.baseline(2)], [a(3),a(4)], 'k--');% baseline end
	plot([options.searchwindow(1),options.searchwindow(1)], [a(3),a(4)], 'r--');% search window start
	plot([options.searchwindow(2),options.searchwindow(2)], [a(3),a(4)], 'r--');% search window end
	if ~isempty(onset)
	    plot([onset,onset], [a(3),a(4)], 'b-');                       % onset detected?
	end
        ylabel('EMG amplitude, mV');
	axis([-options.baseline(1).*1.05,100,a(3),a(4)]);                 % reset axis
	    
	subplot(2,1,2);                                                   % bottom = the diagnostics
	hold on;
	xrange = options.searchwindow(1) : (1000./samplehz) : options.searchwindow(2);% the x-axis time labels
	plot([a(1),a(2)], [base.mean+(options.climit.*base.SD), base.mean+(options.climit.*base.SD)],'k--');% upper detection limit
        plot([a(1),a(2)], [base.mean-(options.climit.*base.SD), base.mean-(options.climit.*base.SD)],'k--');% lower detection limit
	plot(xrange,stats.uppersum, 'b-');                                % the upper cumulative sum
	plot(xrange,stats.lowersum, 'r-');                                % the upper cumulative sum
	a = axis;
	if ~isempty(stats.upper)
	    plot([stats.upper,stats.upper], [0,a(4)], 'b-');              % upper onset detected
	end
	if ~isempty(stats.lower)
	    plot([stats.lower,stats.lower], [0,a(3)], 'r-');              % lower onset detected
        end
	axis([-options.baseline(1).*1.05,100,a(3),a(4)]);                 % reset axis
        xlabel('Time after TMS, ms');
        ylabel('Cumulative sum, mV');
    end
end
