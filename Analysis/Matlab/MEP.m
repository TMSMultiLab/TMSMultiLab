%% measure a single MEP (or mean MEP) from a one [or two] dimensional time series
% [mep]=MEP(data,samplehz,pulsetime[,options])
% inputs required:
% data (nx1) - single channel of time series data; rotated if one row
% samplehz - single number for sampling frequency in samples per second, eg, 2000
% pulsetime - time in ms when electrical or magnetic stimulus was given, relative to the start of data, eg 1000+
% options
%	options.auc - how to calcualte the area under the curve? default: 'threshold'
%		'threshold': use same criteria as threshold (above)
%		'window': use the MEP search window
%	options.average (false,true) default: false
%		false:	treat each input sweep separately, add multiple sweeps into structured output
%		true:	average sweeps before analysis (rectification is applied before averaging)
%	options.baseline [ms,ms] - when is the baseline period, in ms *before* pulsetime? default: all data until pulsetime
%		one value: onset only
%		two values: onset & offset
%		NaN for no limit
%	options.biometrics.age (eg 0:150) - participant age in years
%	options.biometrics.armlength (eg 50:150) - participant arm length in cm
%	options.biometrics.height (eg 100:300) - participant height in cm
%	options.demean (true,false) - remove DC offset (mean of baseline); applied before rectification; default: true
%	options.filter - [low,high] - low and high-pass filters used on the data already
%   options.flip (false,true) - flip the y-values if the minimum peak is earlier than the maximum? default: false
%	options.mirror (true,false) - default: true
%		true:	whatever data are extracted from the post-TMS period (amplitude, AUC, latency), get the same data as a control from the pre-TMS period, with mirrored windows and algorithms
%		false:	only extract post-TMS data
%	options.muscle - label for the muscle recorded ('FDI','APB', etc)
%	options.mvc - level of MVC (0:1). complete rest=0, 10%MVC = 0.1
%	options.plot (false,true) - plot a figure of the anotated MEP? default: false
%	options.rectify (false,true) - rectify data before measuring? default: false
%	options.sweeps - (1:inf) - how many individual trials averaged in the input trace? default: 1
%	options.threshold.type ('abs','p','peak','rel','rms','sd') - what is the MEP (onset/offset) threshold type? default: peak
%		'abs': absolute, in mV (eg 0.05) relative to baseline mean (~zero)
%		'p': p-value threshold relative to baseline mean & standard deviation
%       'peak': max and min within the window (default)
%		'rel': relative to mean rectified
%		'rms': relative to root-mean-squared of baseline
%		'sd': relative to mean & standard deviation of baseline
%	options.threshold.criterion (0:Inf) - numerical criterion for options.threshold, default depends on threshold.type
%		(abs): number of mV - eg 0.05 for 50uV
%		(p): p-value cutoff - eg 0.05 or 0.01
%		(rel): proportion - eg 1.5 for 50% above baseline, 1.0 for equal to baseline
%		(rms): proportion - eg 2.0 for 2x RMS
%		(sd):  number of SDs - eg 2.33
%	options.threshold.direction, default='peak'
%		'peak' - start at the absolute maximum peak then look backwards for onset and forwards for offset
%		'window' - start at the beginning of the MEP window and go forward (default)
%	options.threshold.duration [on,off] (eg 10,0) - for how long in ms does the criterion need to be exceeded to count an onset/offset? default: [0,0] (ie first sample)
%		one value:	use same duration for onset and offset
%		two values:	use first for onset, second for offset
%	options.threshold.proportion (0:1) - what proportion of samples within the specified duration must be above criterion? default: 1 (all of them)
% 	options.window [ms,ms,ms] - when to search for MEP peaks; default: 10-50ms
%		one value: onset only
%		two values: onset & offset
%		three values: onset, offset, & max MEP duration (offset window extended if necessary)
%		NaN for no limit
%
%% OUTPUTS
% mep	structure with many fields; when analysis is mirrored, first value gives post-TMS data, second pre-TMS data
%	mep.amp		peak-to-peak amplitude, mV
%	mep.auc		area under curve, mV.ms
%	mep.minamp	amplitude of lowest peak, mV
%	mep.minlat	latency of lowest peak, ms
%	mep.maxamp	amplitude of highest peak, mV
%	mep.maxlat	latency of highest peak, ms
%	mep.onset	latency of mep onset, ms
%	mep.offset	latency of mep offset, ms
%	mep.duration	duration of mep, from onset to offset, ms
%	mep.baseline	mean, max, min, SD, RMS of baseline
%	mep.flipped	(false,true) - was data reversed during analysis?

% height / arm-length corrections - re-calculate mep latencies by scaling to population average height / arm length

function [mep,options]=MEP(data,samplehz,pulsetime,options)

    %% check input parameters________________________________________________
    if nargin==0
        error('data not given, samplehz not given, pulsetime not given');
    elseif nargin==1
        error('samplehz not given, pulsetime not given');
    elseif nargin==2
        error('pulsetime not given');
    end
    if isempty(samplehz)
        error('samplehz not given');
    end
    if isempty(pulsetime)
        error('pulsetime not given');
    end
    data=squeeze(data);
    if ~ismatrix(data)
        error('data has too many dimensions, should be nx1 array');
    end
    if size(data,1)==1
        data=data';
    end
    if pulsetime>size(data,1)/(samplehz./1000)
        error('pulsetime is later than the end of the data');
    end
    
    %% process the pulsetime_________________________________________________
    endtime=(size(data,1)./(samplehz./1000));
    
    %% load all the default options__________________________________________
    if nargin==3
        options.auc='threshold';                                            % base AUC calculation on the threshold criteria
        options.average=false;                                              % the data are not an average
	    options.baseline=[pulsetime-(1000./samplehz),0];                    % use all available baseline data
	    options.demean=true;                                                % remove the mean across the baseline
        options.flip=false;                                                 % don't invert data if first peak is not the positive max
	    options.mirror=true;                                                % calculate both the post-TMS data and a mirrored portion of data pre-TMS
	    options.rectify=false;                                              % don't rectify the data
	    options.sweeps=1;                                                   % only a single sweep has been passed
	    options.threshold.type='peak';                                      % use the min and max peaks to calculate an MEP (no other options, no onset/offset)
        options.window=[NaN,NaN,NaN];                                       % default window to find MEPs after TMS
        % plotting options
	    options.plot=true;                                                  % plot the data
        options.figure=1;                                                   % Figure 1
        options.subplot=NaN;                                                % no sub-plots
        options.title='';                                                   % no title
    end

    %% check the options.baseline parameter__________________________________
    if isfield(options,'baseline') && ~isempty(options.baseline)
        if abs(options.baseline(1))>pulsetime
	        error('options.baseline starts before the data (or pulsetime is too late)');
        end
	    if isempty(options.baseline(2))
	        options.baseline(2)=0;
	    end
	    if options.baseline(2)<options.baseline(1)
	        error('options.baseline end is before the start');
	    end
    else
        options.baseline=[-pulsetime,0];
    end

    %% check the options.demean parameter____________________________________
    if ~isfield(options,'demean') || ~islogical(options.demean)
        options.demean=true;
    end
    
    %% check the options.mirror parameter____________________________________
    if ~isfield(options,'mirror') || ~islogical(options.mirror)
        options.mirror=true;
    end
    
    %% check the options.plot parameter______________________________________
    if ~isfield(options,'plot') || ~islogical(options.plot)
        options.plot=true;
        if ~isfield(options,'figure') || ~isfinite(options,'figure')
            options.figure=1;
        end
        if ~isfield(options,'subplot')
            options.subplot=NaN;
        end
        if ~isfield(options,'title')
            options.title='';
        end
    end
    
    %% check the options.threshold parameters________________________________
    if ~isfield(options.threshold,'type')
        options.threshold.type='peak';                                      % default = use peak-to-peak and no onset/offset
    else
        switch options.threshold.type
            case 'abs'
                if ~isfield(options.threshold,'criterion')
    	            options.threshold.criterion=0.05;                       % absolute threshold for MEP onset = 0.05mV
                end
            case 'p'
                options.threshold.criterion=0.001;                          % p-value threshold = 0.001 (Z=3.09)
            case 'peak'
            case 'rel'
                options.threshold.criterion=2;                              % relative criterion = 2 times the max / min baseline
            case 'rms'
            case 'sd'
                if ~isfield(options.threshold,'criterion')
    	            options.threshold.criterion=3.09;                       % SD threshold = 3.09 (p<.001)
                end
        end
        if ~isfield(options.threshold,'duration')
	        options.threshold.duration=[0,0];                               % duration that threshold must be exceeded, default = first sample
        end
        if ~isfield(options.threshold,'proportion')
	        options.threshold.proportion=1;                                 % proportion of samples that must be above/below threshold within the above duration, default = all
        end
    end

    %% check the options.window parameter____________________________________
    if isfield(options,'window') && ~isempty(options.window)
        if ~ismatrix(options.window)
            error('options.window has too many dimensions, should be 2x1 array');
        end
        if numel(options.window)>3
            error('options.window has the too many elements, max 3x1 array');
        end
        if isempty(options.window(2))
            options.window(2)=50;
        end
        if isempty(options.window(3))
            options.window(3)=NaN;
        end
	    if ~isfinite(options.window(1))
	        options.window(1)=0; % start immediately after TMS pulse
	    end
	    if ~isfinite(options.window(2))
	        options.window(2)=endtime-pulsetime; % look until the end of the data
	    end
    else
        options.window=[10,50,NaN]; % default window
    end

    %% convert parameters to sample numbers (_s, confusingly)________________
    pulsetime_s=round(pulsetime.*(samplehz./1000));
    endtime_s=round(endtime.*(samplehz./1000));
    options.baseline_s=round((pulsetime+options.baseline).*(samplehz./1000));
    options.baseline_s(1)=options.baseline_s(1)+1;                          % minimum = first sample
    options.window_s=round(options.window.*(samplehz./1000));

    %% MEASURE THE MEP_______________________________________________________
    sweep=nan(size(data,1),options.mirror+1);                               % pre-allocate 
    
    for n=1:options.mirror+1 % first use post-pulse data, then use pre-pulse data as a control
        sweep(:,n)=data;                                                    % use the post-pulse data
        if n~=1
            sweep(pulsetime_s+1:endtime_s,n)=data(1:pulsetime_s);           % copy pre-pulse data into post-pulse data
        end

        % get start and end sample of analysis window
        start=pulsetime_s+options.window_s(1);
        finish=pulsetime_s+options.window_s(2);
	    
	    % subtract the mean baseline from the data
	    if options.demean
	        sweep(:,n)=sweep(:,n)-mean(sweep(options.baseline_s(1):options.baseline_s(2),n));
	    end
	    
	    %% MEASURE THE BASELINE______________________________________________
	    if n==1
            mep.baseline.mean=mean(sweep(options.baseline_s(1):options.baseline_s(2),n));
            mep.baseline.max=max(sweep(options.baseline_s(1):options.baseline_s(2),n));
            mep.baseline.min=min(sweep(options.baseline_s(1):options.baseline_s(2),n));
            mep.baseline.sd=std(sweep(options.baseline_s(1):options.baseline_s(2),n));
            mep.baseline.rms=sqrt(mean(sweep(options.baseline_s(1):options.baseline_s(2),n).^2));
	    end

        %% GET PEAKS AND FLIP THE DATA?______________________________________
        % get min and max peaks and latencies
        [mep.minamp(n),mep.minlat(n)]=min(sweep(start:finish,n));
        [mep.maxamp(n),mep.maxlat(n)]=max(sweep(start:finish,n));
    
        % make the first MEP peak the positive one
        if mep.minlat(n)<mep.maxlat(n)
            sweep(:,n)=-sweep(:,n);
            [mep.minamp(n),mep.minlat(n)]=min(sweep(start:finish,n));
            [mep.maxamp(n),mep.maxlat(n)]=max(sweep(start:finish,n));
            mep.flipped(n)=true;
        else
	        mep.flipped(n)=false;
        end
        
        % FLIPPING IS PROBLEMATIC WHEN ANOTHER LARGE SIGNAL WITHIN WINDOW
	    % DECISION TO FLIP COULD BE BASED ON MEAN MEP DIRECTION PER PARTICIPANT AND MUSCLE?
	    % OR SET TIGHTER WINDOWS BASED, EG, ON PARTICIPANT HEIGHT AND MUSCLE?
        
        %% GET ONSET / OFFSET LATENCIES______________________________________
        switch options.threshold.type
            case 'peak'
            case 'sd'
                supra(:,n)=sweep(:,n)./sweep(:,n);                          % index of ones as long as the sweep, then replace segments with NaN
                supra(abs(sweep(:,n))<(mep.baseline.mean+options.threshold.criterion.*mep.baseline.sd),n)=0;% datapoints outside of the threshold
    
                % which segments are > the threshold.duration?
                if options.threshold.duration>0
                    onsets=find(diff(supra(:,n))==1);                       % find when supra goes from 0 to 1 = onset
                    offsets=find(diff(supra(:,n))==-1);                     % find when supra goes from 1 to 0 = offset
    
                    if onsets(1)>offsets(1)                                 % if first sample is already above threshold
                        onsets=[1;onsets];                                  % add another onset = first sample
                    end
                    if numel(onsets)>numel(offsets)                         % if last sample is still above threshold
                        offsets=[offsets;endtime_s];                        % add another offset = last sample
                    end
                    for o=1:numel(onsets)
                        if (offsets(o)-onsets(o)+1) < (options.threshold.duration.*samplehz/1000) % if segment is too short
                            supra(onsets(o)+1:offsets(o),n)=0;              % remove from consideration
                        end
                    end
                end
                supra(supra(:,n)==0,n)=NaN;                                 % replace zeros with NaN
        end

        %% GET MEP PEAK-TO-PEAK AMPLITUDE____________________________________
        mep.amp(n)=mep.maxamp(n)-mep.minamp(n);
    end
    
    %% convert samples back to ms
    mep.minlat=(mep.minlat+options.window_s(1)-1).*(1000./samplehz);
    mep.maxlat=(mep.maxlat+options.window_s(1)-1).*(1000./samplehz);
    
    %% PLOT THE DATA_______________________________________________________
    if options.plot
        xt=-pulsetime+(1000./samplehz):(1000./samplehz):(endtime-pulsetime);% timescale for x-axis

        % where to plot?
        figure(options.figure);
        if isfinite(options.subplot)
            subplot(options.subplot(1),options.subplot(2),options.subplot(3));% choose the sub-plot
        end
	    hold on;
	    plot(xt,sweep(:,1),'b-');                	                        % plot the data

        if strcmp(options.threshold.type,'sd')                              % plot data outside the threshold = potential response
            plot(xt,sweep(:,1).*supra(:,1),'m');
        end
	    
	    % adjust the axis
	    if options.window(2)>(endtime-pulsetime)./2
	        plotxend=endtime-pulsetime; % plot to end of data
	    else
	        plotxend=options.window(2).*2;% plot to twice the window length
	    end
	    axis([-options.window(2),plotxend,-1,1]);
	    axis 'auto y';
	    
	    % plot max and min peaks
	    plotxwidth=options.window(2)./10;
	    plot([mep.maxlat(1)-plotxwidth,mep.maxlat(1)+plotxwidth],[mep.maxamp(1),mep.maxamp(1)],'b:');
	    plot([mep.minlat(1)-plotxwidth,mep.minlat(1)+plotxwidth],[mep.minamp(1),mep.minamp(1)],'b:');
	    
	    % plot pre-pulse baseline levels
	    xrng=[options.baseline(1),options.baseline(2)];
	    plot(xrng,[mep.baseline.mean,mep.baseline.mean],'k-');
        for s=1:5 % standard deviation lines
	        plot(xrng,[mep.baseline.mean+s.*mep.baseline.sd,mep.baseline.mean+s.*mep.baseline.sd],'--','color',[0,0,0]+0.15.*s);
	        plot(xrng,[mep.baseline.mean-s.*mep.baseline.sd,mep.baseline.mean-s.*mep.baseline.sd],'--','color',[0,0,0]+0.15.*s);
        end
	    plot(xrng,[mep.baseline.max,mep.baseline.max],'k:');
	    plot(xrng,[mep.baseline.min,mep.baseline.min],'k:');
	    
        a=axis;
	
	    % plot x- and y-axis lines
	    plot([-options.window(2),options.window(2).*2],[0,0],'color',[0.6,0.6,0.6]);
	    plot([0,0]+(1000./samplehz),[a(3:4)],'color',[0.6,0.6,0.6]);
	    
	    % plot MEP analysis window
	    plot([options.window(1),options.window(1)],[a(3:4)],'g-');
	    plot([options.window(2),options.window(2)],[a(3:4)],'r-');

	    % label and print
        a=axis;
        if strcmp(options.threshold.type,'sd')
            text(0+((a(2)-a(1))*0.075),a(4)-0.025*sum(abs(a(3:4))),['threshold: ',num2str(options.threshold.criterion),' SD, ',num2str(options.threshold.duration),' ms']);
        end
        title(options.title);
	    xticks(-options.window(2):options.window(2)./5:options.window(2).*2);
	    xlabel('time after pulse, ms');
	    ylabel('EMG signal, mV');
	    print('mep.png','-dpng');
    end
end