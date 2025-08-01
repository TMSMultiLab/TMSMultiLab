%% process an array of epoched EMG data during background contraction with stimulus every epoch
% [output, options, bootstrap] = cutaneomotor(data,samplehz,stimtime [,opts])
% inputs required:
% data - (m x n array, where m = samples per epoch; n = repetitions / epochs)
% samplehs - sampling frequency in hz (samples per second)
% stimtime - time in s at which the stimulus was applied (same on every epoch, measured from the start of the epoch; at least 0.05 s (50ms) is required)
% opts - options for this analysis, with default options indicated, include:
%   opts.artefact = true            remove the stimulus artefact; default = true
%   opts.filter = true              filter the data; default = true. can only be applied if opts.artefact is true
%   opts.padding = 10;              padding in samples at start and end of epoch - to exclude filtering effects
%   opts.criterion = 1.96           criterion for onset/offset of a change from baseline; defaults to 1.96 standard deviations from the baseline mean across samples
%   opts.latency = 20               minimum time in ms after the stimulus that cutaneomotor reflexes could occur - e.g., the H-reflex latency
%   opts.duration = []              minimum duration for a significant epoch, in ms; if specified, then bootstrap is not used
%   opts.bootstrap = true;          bootstrap the criterion for the minimum duration of a significant sequence?
%   opts.iterations = 1000;         iterations in the bootstrap
%   opts.plot = true                plot the data
%   opts.figure = 1                 which figure to plot to
%   opts.subplot = [1,1,1]          which subplot to plot to
%   opts.title = [];                title for the figure
%   opts.verbose = false            print info to screen
%
%   Version 1.0
%   01 August 2025
%   Nick Holmes https://github.com/TheHandLab

function [output, opts, bootstrap] = cutaneomotor(data,samplehz,stimtime,opts)

    %% CHECK THE INPUTS__________________________________________________
    if nargin < 3
        error('output=cutaneomotor(data,samplehz,stimtime [,opts]) - requires at least three input arguments');
    end
    if ndims(data)>2
        error('output=cutaneomotor(data,samplehz,stimtime [,opts]) - data must be a two dimensional array, m x n, with m samples per epoch and n epochs');
    end
    if stimtime < 0.05
        error('output=cutaneomotor(data,samplehz,stimtime [,opts]) - stimtime must be at least 0.05s (50ms), or at least 0.1s (100ms) with artefact removal and filtering');
    end
    if nargin == 3                                                          % load the default options
        opts.artefact=true;                                                 % remove the stimulus artefact using spike_artefact.m
        opts.filter=true;                                                   % filter the EMG data using EMG_filter.m
	opts.padding = 10;                                                  % padding in samples at start and end of epoch - to exclude filtering effects
        opts.criterion = 1.96;                                              % use 1.96 standard deviations (across samples of the baseline) as the criterion for creating sequences
        opts.latency = [20,250];                                            % range of times after stimulus to search for changes in EMG signals
        opts.duration = [];                                                 % minimum duration for a significant sequence, in ms (if bootstrapping not used); empty = use bootstrap
        opts.bootstrap = true;                                              % bootstrap the criterion for the minimum duration
        opts.iterations = 1000;                                             % iterations in the bootstrap
        opts.plot=true;                                                     % plot the data
        opts.figure=1;                                                      % to figure 1
        opts.subplot=[1,1,1];                                               % to subplot 1
	opts.title = [];                                                    % no title
        opts.verbose = false;                                               % don't print info to screen
    else
        if ~isfield(opts,'artefact')
            opts.artefact=true;
        end
        if ~isfield(opts,'filter')
            opts.filter=true;                                               % filter the EMG data using EMG_filter.m
        end
        if ~isfield(opts,'padding')
            opts.padding=10;                                                % padding in samples at start and end of epoch - to exclude filtering effects
        end
        if ~isfield(opts,'criterion')
            opts.criterion = 1.96;                                          % use 1.96 standard deviations (across samples of the baseline) as the criterion
        end
        if ~isfield(opts,'latency')
            opts.latency = [20,250];                                        % range of times after stimulus to search for changes in EMG signals
        end
        if ~isfield(opts,'duration')
            opts.duration = [];                                             % minimum duration for a significant epoch, in ms (ignored if bootstrap used)
	    opts.bootstrap = true;                                          % bootstrap the criterion for the minimum duration                                          
        end
        if ~isfield(opts,'bootstrap')
            opts.bootstrap = true;                                          % bootstrap the criterion for the minimum duration
	end
        if isnumeric(opts.duration) & ~isempty(opts.duration)
	    opts.bootstrap = false;                                         % if a value for the duration is given, turn off the boostrap
        end
        if ~isfield(opts,'iterations')
            opts.iterations = 1000;                                         % iterations in the bootstrap
        end                                             
        if ~isfield(opts,'plot')
            opts.plot=true;                                                 % plot the data
        end
        if ~isfield(opts,'figure')
            opts.figure = 1;                                                % to figure 1
        end
        if ~isfield(opts,'subplot')
           opts.subplot=[1,1,1];                                            % to subplot 1
        end
        if ~isfield(opts,'title')
           opts.title=[];                                                   % no title
        end
        if ~isfield(opts,'verbose')
            opts.verbose = false;                                           % don't print output
        end
    end
    if ~opts.artefact && opts.filter
        warning('output=cutaneomotor(data,samplehz,stimtime [,opts]) - filter option can only be applied with artefact removal');
    end
    if stimtime < 0.1 & opts.filter
        error('output=cutaneomotor(data,samplehz,stimtime [,opts]) - stimtime must be at least 100ms to use artefact removal and filtering');
    end

    %% PROCESS THE INPUTS________________________________________________
    stimtimems= stimtime .* 1000;                                           % convert stimtime in s to stimtime in ms (all times now in ms)
    stimtimes = round(stimtimems .* (samplehz ./ 1000));                    % convert stimtime in ms to stimtime in samples
    if stimtimes > size(data,1)
        error('output=cutaneomotor(data,samplehz,stimtime [,opts]) - stimtime is later than the last sample');
    end
    if ~isempty(opts.duration)
        opts.durations = opts.duration .* (samplehz ./ 1000);               % convert to samples
    else
        opts.durations = [];                                                % default duration length
    end
    samples = size(data,1);                                                 % how many samples in the epoch
    reps = size(data,2);                                                    % how many epochs
    output.cleaned = nan(samples,reps,2);                                   % cleaned data; 3rd dimension, 1: absolute data, baseline subtracted; 2: relative data, divided by baseline
    if opts.filter
        base = stimtimems - 20;                                             % with filtering, baseline goes up to 20ms before the stimulus (at least half of the data before stimtime)
    else
        base = stimtimems - 10;                                             % without filtering, baseline goes up to 10ms before the stimulus (at least half of the data before stimtime)
    end
    bases = round(base .* (samplehz ./ 1000));                              % convert baseline in ms to baseline in samples
    opts.latencys = opts.latency .* (samplehz ./ 1000);                     % convert latency in ms to latency in samples
   
    %% PROCESS THE EPOCHS__________________________________________________
    if opts.verbose
        disp(['Processing ',int2str(reps),' sweeps...']);
    end
    for n = 1:reps
        if mod(n,100)==0 && opts.verbose
            disp([' n=',int2str(n),'...']);                                 % progress on screen
        end
        % remove stimulus artefact & filter the data_______________________
        if opts.artefact
            spike=spike_artefact(data(:,n),samplehz,stimtimems);            % remove spike artefact with default options
            if opts.filter
                output.cleaned(:,n,1)=spike.filtered;                       % the filtered data after spike removal
            else
                output.cleaned(:,n,1)=spike.recovered;                      % the raw data after spike removal
            end
        end
        
        % rectify the data_________________________________________________
        output.cleaned(:,n,1) = abs(output.cleaned(:,n,1));                 % make all data positive -> overall change in EMG activity

        % subtract the baseline before the stimulus________________________
        baseM = nanmean(output.cleaned(opts.padding:bases,n,1),1);          % get mean signal before stimtime
        output.cleaned(:,n,2) = output.cleaned(:,n,1) ./ baseM;             % second output type is divided by the baseline (mean baseline = 1)
        output.cleaned(:,n,1) = output.cleaned(:,n,1) - baseM;              % first output type is with the baseline subtracted (mean baseline = 0)
    end
    
    %% BOOTSTRAP THE MINIMUM LENGTH OF SAMPLES FOR A SIGNIFICANT SEQUENCE__
    if opts.bootstrap
        if opts.verbose
            disp(['Bootstrapping the baselines, ',int2str(opts.iterations),' iterations...']);
        end
        output.sequences=[];                                                % empty variable for the significant sequences
        for n = 1:opts.iterations
            if mod(n,100)==0 && opts.verbose
                disp([' i=',int2str(n),'...']);                             % progress on screen
            end
            bootstrap.reps = ceil(rand(reps,1) .* reps);                    % which repetitions to sample from output.cleaned
            bootstrap.mean = squeeze(nanmean(output.cleaned(opts.padding:bases,bootstrap.reps,2),2));% mean across the resampled repetitions - relative data only
            bootstrap.baseM = nanmean(bootstrap.mean);                      % mean across timepoints (should be zero for the first dimension, and close to 1 for the second
            bootstrap.baseSD = nanstd(bootstrap.mean);                      % SD across timepoints
            bootstrap.critical(1) = bootstrap.baseM - opts.criterion .* bootstrap.baseSD;% lower cut-off for non-significant periods
            bootstrap.critical(2) = bootstrap.baseM + opts.criterion .* bootstrap.baseSD;% upper cut-off for non-signficant periods
            bootstrap.significant = bootstrap.mean < bootstrap.critical(1) | bootstrap.mean > bootstrap.critical(2);% which timepoints are significant?
            bootstrap.onsets = 1+find(diff(bootstrap.significant)==1);      % when does a non-significant timepoint become signifant?
            bootstrap.offsets = 1+find(diff(bootstrap.significant)==-1);    % when does a significant timepoint become non-signifant?
            if bootstrap.offsets(1) < bootstrap.onsets(1)
                bootstrap.onsets=[1;bootstrap.onsets];                      % add an onset before the first timepoint
            end
            if numel(bootstrap.onsets) > numel(bootstrap.offsets)
                bootstrap.offsets = [bootstrap.offsets;numel(bootstrap.significant)+1]; % add a final offset after the last timepoint
            end
            output.sequences = [output.sequences;bootstrap.offsets - bootstrap.onsets];% append the new list of significant sequences
        end
        output.sequences = sortrows(output.sequences);                      % sort the 'significant' sequences by length
        output.criterion = output.sequences(round(numel(output.sequences).*0.99));% Criterion = 99th percentile of randomised sequence lengths
    else
        output.criterion = opts.durations;                                  % use the requested sequence length
    end

    %% AVERAGE THE EPOCHS________________________________________________
    output.mean = squeeze(nanmean(output.cleaned,2));                       % average across epochs

    %% GET MEAN AND SD ACROSS MEAN BASELINE SAMPLES______________________
    output.baseM = nanmean(output.mean(opts.padding:bases,:),1);            % mean across timepoints (should be zero for the first dimension, and close to 1 for the second)
    output.baseSD = nanstd(output.mean(opts.padding:bases,:),0,1);          % SD across timepoints
    critical(:,1) = output.baseM - opts.criterion .* output.baseSD;         % lower cut-off for non-significant periods
    critical(:,2) = output.baseM + opts.criterion .* output.baseSD;         % upper cut-off for non-signficant periods
    
    %% FIND SIGNIFICANT PERIODS AFTER MINIMUM LATENCY____________________
    output.significant = nan(size(output.mean,1),1);                        % for storing significant samples
    n=2; % do this for the relative data only (change later?)
        output.significant = output.mean(:,n) < critical(n,1) | output.mean(:,n) > critical(n,2);% which timepoints are initially above the sequence-creation threshold?
        onsets = 1+find(diff(output.significant)==1);                       % when does a non-significant timepoint become signifant?
        offsets = 1+find(diff(output.significant)==-1);                     % when does a significant timepoint become non-signifant?
        if offsets(1) < onsets(1)                                           % if the first offset is before the first onset
            onsets=[1;onsets];                                              % add an onset before the first timepoint
        end
        if numel(onsets) > numel(offsets)                                   % if there are more onsets than offsets
            offsets = [offsets;numel(output.significant)+1];                % add a final offset after the last timepoint
        end	
        output.significant(1:stimtimes+opts.latencys(1)) = 0;               % remove the too-early periods
	output.significant(stimtimes+opts.latencys(2):end) = 0;             % remove the too-late periods
    
        for m = 1:numel(onsets)                                             % for each remaining significant sequence
	    if offsets(m) - onsets(m) < output.criterion                    % if there are too-few samples in this significant sequence
	        if offsets(m) > numel(output.significant)                   % the last offset may be one sample longer than the data
		    o = numel(output.significant);                          % use the last data sample
		else
		    o = offsets(m);                                         % use the offset sample
		end
	        output.significant(onsets(m):o) = 0;                        % remove this significant sequence
            end
        end
    output.significant=logical(output.significant);                         % convert back to a logical index
    
    %% DESCRIBE THE SIGNIFICANT PERIODS FOUND____________________________
    output.onsets = 1+find(diff(output.significant)==1);                    % when does a non-significant timepoint become signifant?
    output.offsets = 1+find(diff(output.significant)==-1);                  % when does a significant timepoint become non-signifant?
    if output.offsets(1) < output.onsets(1)                                 % if the first offset is before the first onset
        output.onsets=[1;output.onsets];                                    % add an onset before the first timepoint
    end
    if numel(output.onsets) > numel(output.offsets)                         % if there are more onsets than offsets
        output.offsets = [output.offsets;size(output.significant,1)+1];     % add a final offset after the last timepoint
    end
    output.segments = struct('direction',[],'onset',[],'offset',[],'peak',[],'mean',cell(numel(output.onsets),1));% cell array with one per period each containing a structure with elements
    for n = 1:numel(output.onsets)
        if output.mean(output.onsets(n),2) > 1
            output.segments(n).direction = 'increase';                      % increased EMG signal
	else
            output.segments(n).direction = 'decrease';                      % decreased EMG signal
	end
	output.segments(n).onset = (output.onsets(n) - stimtimes) ./ (samplehz ./ 1000);% onset of significant sequence in ms
	output.segments(n).offset = (output.offsets(n) - stimtimes) ./ (samplehz ./ 1000);% offset of significant sequence in ms
	if output.segments(n).direction == 'increase'
            output.segments(n).peak=max(output.mean(output.onsets(n):output.offsets(n),2));% peak change peak change (max for increases, min for decreases)
	else
            output.segments(n).peak=min(output.mean(output.onsets(n):output.offsets(n),2));% peak change peak change (max for increases, min for decreases)
	end
	output.segments(n).mean=mean(output.mean(output.onsets(n):output.offsets(n),2));% average change (mean between onset and offset)
    end
    
    %% PLOT THE DATA_____________________________________________________
    if opts.plot
        xrange = [-stimtimems + (1000 ./ samplehz) : (1000 ./ samplehz) : (samples ./ (samplehz ./ 1000)) - stimtimems]'; % values for the x-axis
        figure(opts.figure);
        subplot(opts.subplot(1),opts.subplot(2),opts.subplot(3));
        hold on;
        plot([0,0],[0,2],'k-');                                             % stimulus time
        plot([xrange(1),xrange(end)],[output.baseM(2),output.baseM(2)],'k-');% baseline
        plot([xrange(1),xrange(end)],[critical(2,1),critical(2,1)],'k:');   % lower cut-off
        plot([xrange(1),xrange(end)],[critical(2,2),critical(2,2)],'k:');   % upper cut-off
        plot(xrange,output.mean(:,2),'b-');                                 % plot the data on top
        sig = output.significant;                                           % which points to highlight
        X=xrange;
        X(output.significant==0) = NaN;                                     % replace with NaN to allow discontinuous plotting
        Y=output.mean(:,2);
        Y(output.significant==0) = NaN;                                     % replace with NaN to allow discontinuous plotting
        plot(X,Y,'r-');                                                     % plot the significant data in red on top
	text(50,0.1,['Criterion: ',int2str(output.criterion),' samples (',num2str(output.criterion ./ (samplehz ./ 1000),3),' ms) > baseline mean ± ',num2str(opts.criterion,3),'SD'],'Color','r');% print the criterion for significance
        xlabel('Time after stimulus, ms');                                  % x axis label
        ylabel('Relative EMG signal, A.U.');                                % y axis label
	title(opts.title);                                                  % add title to the graph
        axis([-50,opts.latency(2),0,2]);                                    % focus on the most-relevant ranges
	set(gcf,'Position',[50,50,1200,800]);                               % set figure size before saving
    end
end