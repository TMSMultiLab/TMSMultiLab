%% FIT A CURVE TO A SET OF DOSE-RESPONSE DATA______________________________
function [fitresult, gof, options] = dose_response(dose, response, options)

    %% CHECK THE INPUTS____________________________________________________
    if nargin < 2
        error('[output, fitresult, options] = dose_response(dose, response [, options]) - at least 2 input arguments required (dose, response)');
    end
    if numel(dose) ~= numel(response) | min(size(dose)) > 1 | min(size(response)) > 1
        error('[output, fitresult, options] = dose_response(dose, response [, options]) - dose and response arguments must be mx1 arrays of the the same size');
    end
    if size(dose,2) > 1
        dose = dose';
    end
    if size(response,2) > 1
        response = response';
    end

    %% SET DEFAULT OPTIONS_________________________________________________
    if nargin == 2 || isempty(options)                                      % if no options passed in
        options.plot = true;                                                % default = plot the data
        options.figure = 1;                                                 % default figure = 1
        options.subplot = [1,1,1];                                          % default subplot = 1 (none)
        options.baseline = 0;                                               % default baseline = 0
        options.baseline_force = false;                                     % force baseline when fit fails; use options.baseline as fixed value
        options.control = [];                                               % default control = empty
        options.xrange = [min(dose),max(dose)];                             % range of data
        options.maxrange = [min(dose) - 0.1.*abs(min(dose)), max(dose) + 0.1.*abs(max(dose))]; % maximum range = 10% either side
        options.reference = NaN;                                            % no default reference point (e.g., a threshold)
        options.xlabel = 'Dose, A.U.';                                      % default x-label
        options.ylabel = 'Response, A.U.';                                  % default y-label
        options.title = 'Dose - response curve';                            % default title
        options.annotate = true;                                            % annotate the graph
        options.fit = 'Sigmoidal';                                          % default (only, at present) curve fit
	options.log = false;                                                % the data should be log-transformed
        options.lowerbounds = [-Inf, -Inf, -Inf, -Inf];                     % minimum possible values for a, b, c, d
        options.upperbounds = [Inf, Inf, Inf, Inf];                         % maximum possible values for a, b, c, d
    end

    %% ADD MISSING DEFAULTS________________________________________________
    if ~isfield(options,'plot')
        options.plot = true;                                                % default = plot the data
    end
    if ~isfield(options,'figure')
        options.figure = 1;                                                 % default figure = 1
    end
    if ~isfield(options,'subplot')
        options.subplot = [1,1,1];                                          % default subplot = 1 (none)
    end
    if ~isfield(options,'baseline')
        options.baseline = 0;                                               % default baseline = 0
    end
    if ~isfield(options,'baseline_force')
        options.baseline_force = false;                                     % default = don't force baseline
    end
    if ~isfield(options,'control')
        options.control = [];                                               % default control = empty
    end 
    if ~isfield(options,'xrange')
        options.xrange = [min(dose),max(dose)];                             % range of data
    end
    if ~isfield(options,'maxrange')
        options.maxrange = [min(dose) - 0.1.*abs(min(dose)), max(dose) + 0.1.*abs(max(dose))]; % maximum range = 10% either side
    end
    if ~isfield(options,'reference')
        options.reference = NaN;                                            % default subplot = 1 (none)
    end  
    if ~isfield(options,'xlabel')
        options.xlabel = 'Dose, A.U.';                                      % default x-label
    end  
    if ~isfield(options,'ylabel')
        options.ylabel = 'Response, A.U.';                                  % default y-label
    end  
    if ~isfield(options,'title')
        options.title = 'Dose - response curve';                            % default title
    end  
    if ~isfield(options,'annotate')
        options.annotate = true;                                            % default = annotate the graph
    end  
    if ~isfield(options,'fit')
        options.fit = 'Sigmoidal';                                          % default (only, at present) curve fit
    end
    if ~isfield(options,'log')
        options.log = false;                                                % default = linear fit
    end 
    if ~isfield(options,'lowerbounds')
        options.lowerbounds = [-Inf, -Inf, -Inf, -Inf];                     % minimum possible values for a (baseline), b (growth), c (midpoint), d (plateau)
    end  
    if ~isfield(options,'upperbounds')
        options.upperbounds = [Inf, Inf, Inf, Inf];                         % maximum possible values for a (baseline), b (growth), c (midpoint), d (plateau)
    end  

    %% LOG AND SUBTRACT BASELINE FROM DATA_________________________________
    if isempty(options.control)                                             % no control values to subtract
        if options.log
            response = log10(response);                                     % log the data
        end
    else                                                                    % subtract control values before fitting  
        if options.log
	
	    %%%%% adjust this to something more sensible - halfway between second-lowest and zero %%%%%
	    options.control(options.control==0) = 0.01;                     
	    
            response = log10(response) - log10(options.control);
	else
	    response = response - options.control;
        end
    end

    %% PLOT THE DATA_______________________________________________________
    if options.plot
        figure(options.figure);                                             % set the figure
        if numel(options.subplot)==3
            subplot(options.subplot(1), options.subplot(2), options.subplot(3));% set the subplot
        elseif numel(options.subplot==1)
            nexttile(options.subplot);
        end
        hold on;
        plot([options.maxrange], [options.baseline,options.baseline], 'k-');% plot baseline
        plot(dose, response,'b*');                                          % plot the data
        xlabel(options.xlabel);                                             % x-axis label
	if options.log
            ylabel(['Log10(',options.ylabel,')']);                          % y-axis label
	else
            ylabel(options.ylabel);                                         % y-axis label
	end
        title(options.title);                                               % graph title
    end

    % subtract the minimum from the curve before fitting___________________
    if min(response)<0
        options.shift = min(response) - ((max(response) - min(response)) ./ 100);% shift by the lowest point on the curve + 1% of the range
	
	% adjust upper and lower bounds depending on the shift
	if options.log
            options.lowerbounds = [0.05, 0.05, 0.05, 0.05];                 % minimum possible values for a, b, c, d
	end
    else
        options.shift = 0;
    end

    %% FIT THE DATA________________________________________________________
    % f(x) = d + (a-d) / (1 + (x/c)^b)                                      fit a 4-parameter logistic model (sigmoidal curve) to M-WAVE data
    % a = initial baseline horizontal (y)
    % b = growth rate (y/x)
    % c = midpoint between a and d = half way up the curve (x)
    % d = final plateau (y)

    % Set up fittype and options:
    ft = fittype( 'd+(a-d)/(1+(x/c)^b)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = options.lowerbounds;
    opts.Upper = options.upperbounds;
    if options.baseline_force                                               % force fit through the baseline
        opts.Lower(1) = options.baseline;
	opts.Upper(1) = options.baseline;
    end
    if options.log                                                          % convert to log scale
        opts.Lower = log10(opts.Lower);
	%opts.Upper = log10(opts.Upper);
    end
    tmp = sortrows([dose, response - options.shift]);                       % sort data in order of increasing dose, and subtract the shift if needed   
    opts.StartPoint = [options.baseline, options.upperbounds(2)./2, mean(options.xrange), max(response)];% min expected EMG noise; guess; 3/4 up the intensity scale; max MEP

    [xData, yData] = prepareCurveData(tmp(:,1), tmp(:,2));                  % not sure what this does...

    % fit model to data____________________________________________________
    [fitresult, gof] = fit(xData, yData, ft, opts);                         % fit the data
    yfit = fitresult.d + (fitresult.a - fitresult.d) ./ (1 + (xData ./ fitresult.c).^fitresult.b );% compute the expected y-values

    % plot the fit_________________________________________________________
    if options.plot
        if options.annotate
            plot([options.maxrange], [fitresult.a, fitresult.a] + options.shift, 'r--');% plot fitted baseline
	    plot([options.maxrange], [fitresult.d, fitresult.d] + options.shift, 'r--');% plot fitted plateau
	    plot([fitresult.c, fitresult.c], [fitresult.a, mean([fitresult.a, fitresult.d])] + options.shift, 'r--');% plot fitted halfway point
            plot([options.maxrange(1), fitresult.c], [mean([fitresult.a, fitresult.d]), mean([fitresult.a, fitresult.d])] + options.shift, 'r--');% plot to the fitted halfway point
        end
        plot(tmp(:,1), yfit + options.shift, 'r-','LineWidth',1.5);         % add fit to the plot
        a=axis;
        offset(1) = abs(max(dose) - min(dose)) .* 0.017;                    % shift text by 1.7% of the x range
        offset(2) = abs(a(4) - a(3)) .* 0.075;                              % shift each line of text by 7.5% of the y range
        if options.annotate
	    if options.log
	       txt = '  log(y) = ';
	    else
	       txt = '  y = ';
	    end
	    txt = [txt, num2str(fitresult.d ,3),' + (',num2str(fitresult.a,3),'-',num2str(fitresult.d,3),') / (1+ ('];
	    if options.shift~=0
	        txt = [txt, '(x ',num2str(options.shift,'%+.3G'),')'];
	    else
	        txt = [txt, 'x '];
	    end
	    txt = [txt, ' / ',num2str(fitresult.c,3),')^{',num2str(fitresult.b,3),'}, r^2 = ',num2str(gof.adjrsquare,3)];
            text(options.maxrange(1)+offset(1),a(4)-(offset(2).*1), txt, 'Color','r');
            text(options.maxrange(1)+offset(1),a(4)-(offset(2).*2), ['  plateau = ',num2str(fitresult.d + options.shift, 3)], 'Color','k','FontSize', 9);% plateau point on Y-axis
            text(options.maxrange(1)+offset(1),a(4)-(offset(2).*3), ['  baseline = ',num2str(fitresult.a + options.shift, 3)], 'Color','k','FontSize', 9);% baseline point on Y-axis
            text(options.maxrange(1)+offset(1),a(4)-(offset(2).*4), ['  growth rate = ',num2str(fitresult.b, 3)], 'Color','k','FontSize', 9);% baseline point on Y-axis
            text(options.maxrange(1)+offset(1),a(4)-(offset(2).*5), ['  50% response = (',num2str(fitresult.c, 3),', ',num2str(mean([fitresult.a, fitresult.d]) + options.shift, 3),')'], 'Color','k','FontSize', 9);% 50% points on X and Y axes
        end
        axis([options.maxrange, a(3), a(4)]);                               % set to max range
        if isfinite(options.reference)
            plot([options.reference, options.reference], [a(3), a(4)], 'k-');% reference point - e.g., a threshold
        end
    end
end