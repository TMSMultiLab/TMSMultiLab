%% measure & remove the stimulus artefact from a single MEP (or mean MEP)
% [spike,fitresult,gof]=spike_artefact(data,samplehz,stimtime [,artefactwindow] [,recovery] [,scale] [,exclude] [,plot_artefact])
% inputs required:
% data (nx1 array)
% sample frequency in hz (samples per second)
% stimulation time (in ms)
% artefact window [ms,ms], if omitted, uses 0-2ms; this is the non-recoverable spike, and is removed completely
% recovery time (ms), if omitted, uses 20ms
% scale (proportion), multiply the post-artefact fit by this number; if omitted = 1
% exclude [ms,ms], exclude this window (e.g., for an M-wave) from the artefact fit; if omitted = [NaN,NaN];
% uses: EMG_filter
% references:
% SARGE: Erez Y, et al. J Neurosci Methods. 2010. PMID: 20542059 Generalized framework for stimulus artifact removal	
% Fisher et al 2012 - TMS artefact removed by double exponential curve, V(t)=A1.e^[-tk1/T1] - A2.e^[-tk2/T2)
	
function [spike,fitresult,gof]=spike_artefact(data,samplehz,stimtime,artefactwindow,recovery,scale,exclude,plot_artefact)

    %% CHECK INPUT ARGUMENTS_______________________________________________
    if nargin == 1
        error('samplehz not given, tms time (in ms) not given');
    elseif nargin == 2
        error('tms time (in ms) not given');
    end
    if isempty(samplehz)
        error('samplehz not given');
    end
    if isempty(stimtime)
        error('stimtime (in ms) not given');
    end
    if nargin < 4
        artefactwindow = [0,2];                                             % default window to find artefacts after TMS = 0 to 2ms
    end
    if nargin < 5
        recovery = 50;                                                      % default artefact recovery time = 50ms (after MEP latency, so check!)
    end
    if nargin < 6
        scale = 1;                                                          % default = don't rescale the yfit data (e.g., fudge to remove bias caused by M-waves)
    end
    if nargin < 7
        exclude = [NaN,NaN];                                                % default = don't exclude data (e.g., around the M-wave)
    end    
    if nargin < 8
        plot_artefact = 0;                                                  % default = don't plot data
    end
    if plot_artefact ~= 1 && plot_artefact ~= 0
        plot_artefact = 0;                                                  % default = don't plot data
    end
    data = squeeze(data);                                                   % remove singleton dimensions
    if ndims(data) > 2
        error('data has too many dimensions, should be nx1 array');
    end
    if size(data,1) == 1
        data = data';
        %warning('data rotated so that samples are along the first dimension');
    end
    if ndims(artefactwindow) > 2
        error('artefactwindow has too many dimensions, should be 2x1 array');
    end
    if size(artefactwindow,1) == 1
        artefactwindow = artefactwindow';
        %warning('artefactwindow rotated to a 2x1 array');
    end
    if numel(artefactwindow) ~= 2
        error('artefactwindow has too many elements, should be 2x1 array');
    end
    if recovery <= artefactwindow(2)
        error('recovery time is shorter than or equal to the artefact window');
    end
    if numel(exclude)>2
        error('exclude window has too many elements, should be 2x1 array');
    end
    if exclude(1)>exclude(2)
        exclude=exclude([2,1]);
    end
    if sum(isnan(exclude))==1
        error('exclude window contains only one number, should be a 2x1 array of NaNs or numbers');
    end
    if exclude(1)<artefactwindow(2)
        error('exclude window onset is before the end of the artefact window, should be at least 2 samples after the artefact window');
    end
    if exclude(2)>recovery
        error('exclude window offset is after the end of the recovery window, should be at least 2 samples less than the recovery window');
    end

    %% CONVERT INPUTS______________________________________________________
    artefactwindows = round(artefactwindow .* (samplehz./1000));            % convert to samples
    stimtimes = round(stimtime .* (samplehz ./ 1000));                      % convert to samples
    recoverys = round(recovery .* (samplehz ./ 1000));                      % convert to samples
    excludes = round(exclude .* (samplehz ./ 1000));                        % convert to samples
    samples = size(data,1);                                                 % total number of samples
    duration = 1000 .* samples ./ samplehz;                                 % sampling duration in ms
    xrange = -stimtime : (1000 ./ samplehz) : duration - stimtime - (1 ./ samplehz);% x-axis range for the data

    %% REMOVE DC OFFSET (mean of time up to artefact window)_______________
    data=data-mean(data(1:stimtimes-1));
    
    %% SAVE ARTEFACT SPIKE in raw data, within Xms of stimulus onset_______
    spike.full = zeros(size(data,1),1);                                    % for the full spike model
    spike.full(stimtimes + artefactwindows(1) : stimtimes + artefactwindows(2)) = data(stimtimes + artefactwindows(1) : stimtimes + artefactwindows(2));

    %% REMOVE NON-RECOVERABLE SPIKE FROM THE DATA__________________________
    clean=data;
    clean(stimtimes + artefactwindows(1) : stimtimes + artefactwindows(2)) = NaN;

    %% GET DIFFERENCE BETWEEN DATA AND SPIKE SO FAR = RESIDUAL ARTEFACT____
    spike.recovery = data - spike.full;
    
    %% EXCLUDE SEGMENT FROM THE FIT?_______________________________________
    if sum(isnan(exclude))==0                                              % if there is a segment to exclude
        start = stimtimes + excludes(1) + 1;                               % start of interpolation window
        finish = stimtimes + excludes(2) - 1;                              % end of interpolation window
	starty = mean(spike.recovery(start-2:start+2));                    % start point for interpolation
	finishy = mean(spike.recovery(finish-2:finish+2));                 % finish point for interpolation
	delta = (finishy - starty) ./ (diff(excludes)-2);                  % mean change between samples within the window
	if size(starty:delta:finishy,2) < size(start:finish,2)
	    tmp=[starty : delta : finishy,finishy];
	else
	    tmp=starty : delta : finishy;
	end
	if delta~=0
            spike.recovery(start:finish) = tmp;                            % interpolate the data
	else
	    spike.recovery(start:finish) = spike.recovery(start);
	end
    end	

    %% FIT DOUBLE EXPONENTIAL CURVE TO THE RESIDUAL ARTEFACT_______________
    % f(x) = a*exp(b*x) + c*exp(d*x)
    x = xrange(stimtimes + artefactwindows(2) +1 : stimtimes + recoverys)';% the x data potentially affected by the artefact
    y = spike.recovery(stimtimes + artefactwindows(2) +1 : stimtimes + recoverys);% the EMG data affected by the artefact   
    [ xData, yData ] = prepareCurveData( x, y );                            % not sure what this does
   
    %% SET UP FIT TYPE AND OPTIONS_________________________________________
    ft = fittype( 'exp2' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0  0.1  0  0];

    %% FIT MODEL TO DATA___________________________________________________
    [fitresult, gof] = fit( xData, yData, ft, opts );                       % fit the model
    spike.yfit = fitresult.a .* exp(fitresult.b .* x) + fitresult.c .* exp(fitresult.d .* x); % get fitted y values

    %% LAST POINT IN YFIT SHOULD BE ZERO TO AVOID STEP-ARTEFACTS___________
    spike.envelope=sqrt((max(xData)-xData) ./ max(xData));                  % multiply model by this to diminish effect of later timepoints (reduces data removal)
    spike.yfit = spike.yfit .* spike.envelope;

    %% ADD FITTED DATA TO THE SPIKE RECONSTRUCTION_________________________
    spike.full(stimtimes + artefactwindows(2) +1 : stimtimes + recoverys) = spike.yfit .* scale; % multiply by scaling factor

    %% REMOVE SPIKE ARTEFACT FROM DATA_____________________________________
    spike.recovered = data - spike.full;

    %% INTERPOLATE BETWEEN THE LAST GOOD SAMPLE AND FIRST GOOD RECONSTRUCTED SAMPLE
    start = stimtimes + artefactwindows(1) - 1;
    finish = stimtimes + artefactwindows(2) + 1;
    delta = (spike.recovered(finish) - spike.recovered(start)) ./ (diff(artefactwindows) + 2);
    if delta~=0
        spike.recovered(start:finish) = spike.recovered(start) : delta : spike.recovered(finish);% interpolate the data
    else
        spike.recovered(start:finish) = spike.recovered(start);
    end

    %% BANDPASS FILTER THE DATA_____________________________________________
    padded = [repmat(spike.recovered(1),0.1.*samplehz,1);spike.recovered;repmat(spike.recovered(end),0.1.*samplehz,1)];% pad with 100ms to avoid filter roll-off
    filtered=EMG_filter(padded,samplehz,1,20,500);                          % filter 20-500Hz
    spike.filtered=filtered(0.1.*samplehz+1:end-0.1.*samplehz);             % remove the pad
	
    %% PLOT THE DATA________________________________________________________
    if plot_artefact

        % TOP SUBPLOT = RAW DATA____________________________________________
        figure;
        subplot(5,2,1);                                                     % full scale
        hold on;
        plot(xrange,data,'k-');
        ylabel('Raw, mV');
        title('Full epoch');

        subplot(5,2,2);                                                     % zoomed-in
        hold on;
        plot(xrange,data,'k-');
        a=axis;
        axis([artefactwindow(1)-10,recovery+10,a(3),a(4)]);
        title('Zoomed in');

        % NEXT SUBPLOT = SPIKE DATA_________________________________________
        subplot(5,2,3);                                                     % full scale
        hold on;
        plot(xrange,spike.full,'r-');
        ylabel('Spike');

        subplot(5,2,4);                                                     % zoomed-in
        hold on;
        plot(xrange,spike.full,'r-');
        a=axis;
        axis([artefactwindow(1)-10,recovery+10,a(3),a(4)]);
        
        % NEXT SUBPLOT = RAW - SPIKE DATA = RESIDUAL ARTEFACT + RESPONSE TO MODEL
        subplot(5,2,5);                                                     % full scale
        hold on;        
        plot(xrange,spike.recovery,'b-');
        plot(xrange(stimtimes + artefactwindows(2) +1 : stimtimes + recoverys),spike.yfit .* scale,'r-');
        ylabel('Residual');

        subplot(5,2,6);                                                     % zoomed-in
        hold on;
        plot(xrange,spike.recovery,'b-');
        plot(xrange(stimtimes + artefactwindows(2) +1 : stimtimes + recoverys),spike.yfit .* scale,'r-');
        a=axis;
        plot([artefactwindow(1),artefactwindow(1)],[a(3),a(4)],'r-');       % artefact window start
        plot([artefactwindow(2),artefactwindow(2)],[a(3),a(4)],'r-');       % artefact window end
        axis([artefactwindow(1)-10,recovery+10,a(3),a(4)]);                 % rescale axis
        
        % NEXT SUBPLOT = MODELLED FIT AND REMOVAL___________________________
        subplot(5,2,7);                                                     % full scale
        hold on;
        plot(xrange,spike.recovered,'b-');

        ylabel('Removed');
        subplot(5,2,8);                                                     % zoomed-in
        hold on;
        plot(xrange,spike.recovered,'b-');
        a=axis;
        plot([0,0],[a(3),a(4)],'k-');
        axis([artefactwindow(1)-10,recovery+10,a(3),a(4)]);                 % rescale axis

        % LAST SUBPLOT = FILTERED CLEANED DATA______________________________
        subplot(5,2,9);                                                     % full scale
        hold on;
        plot(xrange,spike.filtered,'b-');

        ylabel('Filtered');
        subplot(5,2,10);                                                    % zoomed-in
        hold on;
        plot(xrange,spike.filtered,'b-');
        a=axis;
        plot([0,0],[a(3),a(4)],'k-');
        axis([artefactwindow(1)-10,recovery+10,a(3),a(4)]);                 % rescale axis
        xlabel('Time after stimulus, ms');
    end
end