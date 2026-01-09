%% process an array of epoched EMG data during background contraction with stimulus every epoch
% [output, options, bootstrap] = cutaneomotor(data,samplehz,stimtime [,options])
% inputs required:
% data - (m x n array, where m = samples per epoch; n = repetitions / epochs)
% samplehs - sampling frequency in hz (samples per second)
% stimtime - time in s at which the stimulus was applied (same on every epoch, measured from the start of the epoch; at least 0.05 s (50ms) is required)
% options - options for this analysis, with default options indicated, include:
%   options.artefact = true            remove the stimulus artefact; default = true
%   options.filter = true              filter the data; default = true. can only be applied if options.artefact is true
%   options.padding = 10;              padding in samples at start and end of epoch - to exclude filtering effects
%   options.criterion = 1.96           criterion for onset/offset of a change from baseline; defaults to 1.96 standard deviations from the baseline mean across samples
%   options.significance = 0.99        value for the pth percentile of bootstrapped sequences
%   options.latency = [20,250]         minimum time in ms after the stimulus that cutaneomotor reflexes could occur - e.g., the H-reflex latency
%   options.duration = []              minimum duration for a significant epoch, in ms; if specified, then bootstrap is not used
%   options.bootstrap = true;          bootstrap the criterion for the minimum duration of a significant sequence?
%   options.iterations = 1000;         iterations in the bootstrap
%   options.plot = true                plot the data
%   options.figure = 1                 which figure to plot to
%   options.subplot = [1,1,1]          which subplot to plot to (or tile if a single number)
%   options.title = [];                title for the figure
%   options.annotate = true            annotate the graph
%   options.verbose = false            print info to screen
%
%   Version 1.0
%   01 August 2025
%   Nick Holmes https://github.com/TheHandLab

function [output, options, bootstrap] = cutaneomotor(data, samplehz, stimtime, options)

    %% CHECK THE INPUTS__________________________________________________
    if nargin < 3
        error('output=cutaneomotor(data,samplehz,stimtime [,options]) - requires at least three input arguments');
    end
    if ndims(data)>2
        error('output=cutaneomotor(data,samplehz,stimtime [,options]) - data must be a two dimensional array, m x n, with m samples per epoch and n epochs');
    end
    if stimtime < 0.05
        error('output=cutaneomotor(data,samplehz,stimtime [,options]) - stimtime must be at least 0.05s (50ms), or at least 0.1s (100ms) with artefact removal and filtering');
    end
    if nargin == 3                                                          % load the default options
        options.artefact = true;                                            % remove the stimulus artefact using spike_artefact.m
        options.filter = true;                                              % filter the EMG data using EMG_filter.m
        options.padding = 10;                                               % padding in samples at start and end of epoch - to exclude filtering effects
        options.criterion = 1.96;                                           % use 1.96 standard deviations (across samples of the baseline) as the criterion for creating sequences
	options.significance = 0.99                                         % pth percentile of bootstrapped sequences to accept as 'significant'
        options.latency = [20,250];                                         % range of times after stimulus to search for changes in EMG signals    
        options.duration = [];                                              % minimum duration for a significant sequence, in ms (if bootstrapping not used); empty = use bootstrap
        options.bootstrap = true;                                           % bootstrap the criterion for the minimum duration
        options.iterations = 1000;                                          % iterations in the bootstrap
        options.plot = true;                                                % plot the data
        options.figure = 1;                                                 % to figure 1
        options.subplot = [1,1,1];                                          % to subplot 1
        options.title = [];                                                 % no title
        options.annotate = true;                                            % annotate the graph
        options.verbose = false;                                            % don't print info to screen
    else
        if ~isfield(options,'artefact')
            options.artefact = true;
        end
        if ~isfield(options,'filter')
            options.filter = true;                                          % filter the EMG data using EMG_filter.m
        end
        if ~isfield(options,'padding')
            options.padding = 10;                                           % padding in samples at start and end of epoch - to exclude filtering effects
        end
        if ~isfield(options,'criterion')
            options.criterion = 1.96;                                       % use 1.96 standard deviations (across samples of the baseline) as the criterion
        end
        if ~isfield(options,'significance')
            options.significance = 0.99;                                    % use 99th percentile as cut-off for significant sequences
        end
        if ~isfield(options,'latency')
            options.latency = [20,250];                                     % range of times after stimulus to search for changes in EMG signals
        end
        if ~isfield(options,'duration')
            options.duration = [];                                          % minimum duration for a significant epoch, in ms (ignored if bootstrap used)
            options.bootstrap = true;                                       % bootstrap the criterion for the minimum duration                                          
        end
        if ~isfield(options,'bootstrap')
            options.bootstrap = true;                                       % bootstrap the criterion for the minimum duration
        end
        if isnumeric(options.duration) & ~isempty(options.duration)
        options.bootstrap = false;                                          % if a value for the duration is given, turn off the boostrap
        end
        if ~isfield(options,'iterations')
            options.iterations = 1000;                                      % iterations in the bootstrap
        end                                             
        if ~isfield(options,'plot')
            options.plot = true;                                            % plot the data
        end
        if ~isfield(options,'figure')
            options.figure = 1;                                             % to figure 1
        end
        if ~isfield(options,'subplot')
           options.subplot = [1,1,1];                                       % to subplot 1
        end
        if ~isfield(options,'title')
           options.title = [];                                              % no title
        end
        if ~isfield(options,'annotate')
           options.annotate = true;                                         % annotate the graph
        end
        if ~isfield(options,'verbose')
            options.verbose = false;                                        % don't print output
        end
    end
    if ~options.artefact && options.filter
        warning('output=cutaneomotor(data,samplehz,stimtime [,options]) - filter option can only be applied with artefact removal');
    end
    if stimtime < 0.1 & options.filter
        error('output=cutaneomotor(data,samplehz,stimtime [,options]) - stimtime must be at least 100ms to use artefact removal and filtering');
    end

    %% PROCESS THE INPUTS________________________________________________
    stimtimems = stimtime .* 1000;                                          % convert stimtime in s to stimtime in ms (all times now in ms)
    stimtimes = round(stimtimems .* (samplehz ./ 1000));                    % convert stimtime in ms to stimtime in samples
    if stimtimes > size(data,1)
        error('output=cutaneomotor(data,samplehz,stimtime [,options]) - stimtime is later than the last sample');
    end
    if ~isempty(options.duration)
        options.durations = options.duration .* (samplehz ./ 1000);         % convert to samples
    else
        options.durations = [];                                             % default duration length
    end
    samples = size(data,1);                                                 % how many samples in the epoch
    reps = size(data,2);                                                    % how many epochs
    output.cleaned = nan(samples,reps,2);                                   % cleaned data; 3rd dimension, 1: absolute data, baseline subtracted; 2: relative data, divided by baseline
    if options.filter
        base = stimtimems - 20;                                             % with filtering, baseline goes up to 20ms before the stimulus (at least half of the data before stimtime)
    else
        base = stimtimems - 10;                                             % without filtering, baseline goes up to 10ms before the stimulus (at least half of the data before stimtime)
    end
    bases = round(base .* (samplehz ./ 1000));                              % convert baseline in ms to baseline in samples
    options.latencys = options.latency .* (samplehz ./ 1000);               % convert latency in ms to latency in samples
   
    %% PROCESS THE EPOCHS________________________________________________
    if options.verbose
        disp(['Processing ',int2str(reps),' sweeps...']);
    end
    for n = 1:reps
        if mod(n,100) == 0 && options.verbose
            disp([' n=',int2str(n),'...']);                                 % progress on screen
        end
        % remove stimulus artefact & filter the data_______________________
        if options.artefact
            spike = spike_artefact(data(:,n),samplehz,stimtime);            % remove spike artefact with default options
            if options.filter
                output.cleaned(:,n,1) = spike.filtered;                     % the filtered data after spike removal
            else
                output.cleaned(:,n,1) = spike.recovered;                    % the raw data after spike removal
            end
        end
        
        % rectify the data_________________________________________________
        output.cleaned(:,n,1) = abs(output.cleaned(:,n,1));                 % make all data positive -> overall change in EMG activity

        % subtract the baseline before the stimulus________________________
        baseM = nanmean(output.cleaned(options.padding:bases,n,1),1);       % get mean signal before stimtime
        output.cleaned(:,n,2) = output.cleaned(:,n,1) ./ baseM;             % second output type is divided by the baseline (mean baseline = 1)
        output.cleaned(:,n,1) = output.cleaned(:,n,1) - baseM;              % first output type is with the baseline subtracted (mean baseline = 0)
    end
    
    %% BOOTSTRAP THE MINIMUM LENGTH OF SAMPLES FOR A SIGNIFICANT SEQUENCE__
    if options.bootstrap
        if options.verbose
            disp(['Bootstrapping the baselines, ',int2str(options.iterations),' iterations...']);
        end
        output.sequences = [];                                              % empty variable for the significant sequences
        for n = 1:options.iterations
            if mod(n,100) == 0 && options.verbose
                disp([' i=',int2str(n),'...']);                             % progress on screen
            end
            bootstrap.reps = ceil(rand(reps,1) .* reps);                    % which repetitions to sample from output.cleaned
            bootstrap.mean = squeeze(nanmean(output.cleaned(options.padding:bases,bootstrap.reps,2),2));% mean across the resampled repetitions - relative data only
            bootstrap.baseM = nanmean(bootstrap.mean);                      % mean across timepoints (should be zero for the first dimension, and close to 1 for the second
            bootstrap.baseSD = nanstd(bootstrap.mean);                      % SD across timepoints
            bootstrap.critical(1) = bootstrap.baseM - options.criterion .* bootstrap.baseSD;% lower cut-off for non-significant periods
            bootstrap.critical(2) = bootstrap.baseM + options.criterion .* bootstrap.baseSD;% upper cut-off for non-signficant periods
            bootstrap.significant = bootstrap.mean < bootstrap.critical(1) | bootstrap.mean > bootstrap.critical(2);% which timepoints are significant?
            bootstrap.onsets = 1+find(diff(bootstrap.significant)==1);      % when does a non-significant timepoint become signifant?
            bootstrap.offsets = 1+find(diff(bootstrap.significant)==-1);    % when does a significant timepoint become non-signifant?
	    
            if ~isempty(bootstrap.offsets) & ~isempty(bootstrap.onsets) & bootstrap.offsets(1) < bootstrap.onsets(1)
                bootstrap.onsets = [1;bootstrap.onsets];                    % add an onset before the first timepoint
            end
            if numel(bootstrap.onsets) > numel(bootstrap.offsets)
                bootstrap.offsets = [bootstrap.offsets;numel(bootstrap.significant)+1]; % add a final offset after the last timepoint
            end
            output.sequences = [output.sequences;bootstrap.offsets - bootstrap.onsets];% append the new list of significant sequences
        end
	if ~isempty(output.sequences)
            output.sequences = sortrows(output.sequences);                  % sort the 'significant' sequences by length
            output.criterion = output.sequences(round(numel(output.sequences).*options.significance));% Criterion = pth percentile of randomised sequence lengths
	else
	    output.criterion = size(data,1) - samplehz.*stimtime;           % set criterion to max possible post-stimulus
	end
    else
        output.criterion = options.durations;                               % use the requested sequence length
    end

    %% AVERAGE THE EPOCHS________________________________________________
    output.mean = squeeze(nanmean(output.cleaned,2));                       % average across epochs

    %% GET MEAN AND SD ACROSS MEAN BASELINE SAMPLES______________________
    output.baseM = nanmean(output.mean(options.padding:bases,:),1);         % mean across timepoints (should be zero for the first dimension, and close to 1 for the second)
    output.baseSD = nanstd(output.mean(options.padding:bases,:),0,1);       % SD across timepoints
    critical(:,1) = output.baseM - options.criterion .* output.baseSD;      % lower cut-off for non-significant periods
    critical(:,2) = output.baseM + options.criterion .* output.baseSD;      % upper cut-off for non-signficant periods
    
    %% FIND SIGNIFICANT PERIODS AFTER MINIMUM LATENCY____________________
    output.significant = nan(size(output.mean,1),1);                        % for storing significant samples
    n=2; % do this for the relative data only (change later?)
    output.significant = output.mean(:,n) < critical(n,1) | output.mean(:,n) > critical(n,2);% which timepoints are initially above the sequence-creation threshold?
    onsets = 1+find(diff(output.significant) == 1);                         % when does a non-significant timepoint become signifant?
    offsets = 1+find(diff(output.significant) == -1);                       % when does a significant timepoint become non-signifant?
    if offsets(1) < onsets(1)                                               % if the first offset is before the first onset
        onsets = [1;onsets];                                                % add an onset before the first timepoint
    end
    if numel(onsets) > numel(offsets)                                       % if there are more onsets than offsets
        offsets = [offsets;numel(output.significant)+1];                    % add a final offset after the last timepoint
    end    
    output.significant(1:stimtimes+options.latencys(1)) = 0;                % remove the too-early periods
    output.significant(stimtimes+options.latencys(2):end) = 0;              % remove the too-late periods
    
    for m = 1:numel(onsets)                                                 % for each remaining significant sequence
        if offsets(m) - onsets(m) < output.criterion                        % if there are too-few samples in this significant sequence
            if offsets(m) > numel(output.significant)                       % the last offset may be one sample longer than the data
                o = numel(output.significant);                              % use the last data sample
            else
                o = offsets(m);                                             % use the offset sample
            end
            output.significant(onsets(m):o) = 0;                            % remove this significant sequence
        end
    end
    output.significant = logical(output.significant);                       % convert back to a logical index
    
    %% DESCRIBE THE SIGNIFICANT PERIODS FOUND____________________________
    output.onsets = 1+find(diff(output.significant)==1);                    % when does a non-significant timepoint become significant?
    output.offsets = 1+find(diff(output.significant)==-1);                  % when does a significant timepoint become non-significant?
    if numel(output.onsets)>0
        if output.offsets(1) < output.onsets(1)                             % if the first offset is before the first onset
            output.onsets = [1;output.onsets];                              % add an onset before the first timepoint
        end
        if numel(output.onsets) > numel(output.offsets)                     % if there are more onsets than offsets
            output.offsets = [output.offsets;size(output.significant,1)+1]; % add a final offset after the last timepoint
        end
        output.segments = struct('direction',[],'onset',[],'offset',[],'peak',[],'mean',cell(numel(output.onsets),1));% cell array with one per period each containing a structure with elements
        for n = 1:numel(output.onsets)
            if output.mean(output.onsets(n),2) > 1
                output.segments(n).direction = 'increase';                  % increased EMG signal
            else
                output.segments(n).direction = 'decrease';                  % decreased EMG signal
            end
            output.segments(n).onset = (output.onsets(n) - stimtimes) ./ (samplehz ./ 1000);% onset of significant sequence in ms
            output.segments(n).offset = (output.offsets(n) - stimtimes) ./ (samplehz ./ 1000);% offset of significant sequence in ms
            if output.segments(n).direction == 'increase'
                output.segments(n).peak = max(output.mean(output.onsets(n):output.offsets(n),2));% peak change peak change (max for increases, min for decreases)
            else
                output.segments(n).peak = min(output.mean(output.onsets(n):output.offsets(n),2));% peak change peak change (max for increases, min for decreases)
            end
            output.segments(n).mean = mean(output.mean(output.onsets(n):output.offsets(n),2));% average change (mean between onset and offset)
        end
    end
    
    %% PLOT THE DATA_____________________________________________________
    if options.plot
        xrange = [-stimtimems + (1000 ./ samplehz) : (1000 ./ samplehz) : (samples ./ (samplehz ./ 1000)) - stimtimems]'; % values for the x-axis
        figure(options.figure);
        if numel(options.subplot) == 3
            subplot(options.subplot(1),options.subplot(2),options.subplot(3));
        elseif numel(options.subplot) == 1
            nexttile(options.subplot);                                      % for tiled layout graphs
        end
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
        if options.annotate
            text(50,0.1,['Criterion: ',int2str(output.criterion),' samples (',num2str(output.criterion ./ (samplehz ./ 1000),3),' ms) > baseline mean ± ',num2str(options.criterion,3),'SD'],'Color','r');% print the criterion for significance
        end
        xlabel('Time after stimulus, ms');                                  % x axis label
        ylabel('Relative EMG signal, A.U.');                                % y axis label
        title(options.title);                                               % add title to the graph
        axis([-50,options.latency(2),0,2]);                                 % focus on the most-relevant ranges
    end
end