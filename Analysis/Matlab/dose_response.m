%% FIT A CURVE TO A SET OF DOSE-RESPONSE DATA____________________________
function [fitresult, gof, options] = dose_response(dose, response, options)

    %% CHECK THE INPUTS__________________________________________________
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

    %% SET DEFAULT OPTIONS_______________________________________________
    if nargin == 2 || isempty(options)                                      % if no options passed in
        options.plot = true;                                                % default = plot the data
        options.figure = 1;                                                 % default figure = 1
        options.subplot = [1,1,1];                                          % default subplot = 1 (none)
        options.baseline = 0;                                               % default baseline = 0
        options.xrange = [min(dose),max(dose)];                             % range of data
        options.maxrange = [min(dose) - 0.1.*abs(min(dose)),max(dose) + 0.1.*abs(max(dose))]; % maximum range = 10% either side
        options.reference = NaN;                                            % no default reference point (e.g., a threshold)
        options.xlabel = 'Dose, A.U.';                                      % default x-label
        options.ylabel = 'Response, A.U.';                                  % default y-label
        options.title = 'Dose - response curve';                            % default title
        options.fit = 'Sigmoidal';                                          % default (only, at present) curve fit
        options.lowerbounds = [-Inf,-Inf,-Inf,-Inf];                        % minimum possible values for a, b, c, d
        options.upperbounds = [Inf,Inf,Inf,Inf];                            % maximum possible values for a, b, c, d
    end

    %% ADD MISSING DEFAULTS______________________________________________
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
        options.baseline = 0;                                               % default subplot = 1 (none)
    end    
    if ~isfield(options,'xrange')
        options.xrange = [min(dose),max(dose)];                             % range of data
    end
    if ~isfield(options,'maxrange')
        options.maxrange = [min(dose) - 0.1.*abs(min(dose)),max(dose) + 0.1.*abs(max(dose))]; % maximum range = 10% either side
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
    if ~isfield(options,'fit')
        options.fit = 'Sigmoidal';                                          % default (only, at present) curve fit
    end  
    if ~isfield(options,'lowerbounds')
        options.lowerbounds = [-Inf,-Inf,-Inf,-Inf];                        % minimum possible values for a, b, c, d
    end  
    if ~isfield(options,'upperbounds')
        options.upperbounds = [Inf,Inf,Inf,Inf];                            % maximum possible values for a, b, c, d
    end  

    %% PLOT THE DATA_____________________________________________________
    if options.plot
        figure(options.figure);                                             % set the figure
        subplot(options.subplot(1),options.subplot(2),options.subplot(3));  % set the subplot
        hold on;
        plot([options.maxrange],[options.baseline,options.baseline],'k-');  % plot baseline
        plot(dose,response,'b*');                                           % plot the data
        a=axis;                                                             % get current axis limits
        axis([options.maxrange,a(3),a(4)]);                                 % set to max range
        if isfinite(options.reference)
            plot([options.reference,options.reference],[a(3),a(4)],'k-');   % reference point
        end
        xlabel(options.xlabel);                                             % x-axis label
        ylabel(options.ylabel);                                             % y-axis label
        title(options.title);                                               % graph title
    end

    %% FIT THE DATA______________________________________________________
    % f(x) = d + (a-d) / (1 + (x/c)^b)                                      fit a 4-parameter logistic model (sigmoidal curve) to M-WAVE data
    % a = initial baseline horizontal
    % b = growth rate
    % c = midpoint between a and d = half way up the curve
    % d = final plateau

    % Set up fittype and options:
    ft = fittype( 'd+(a-d)/(1+(x/c)^b)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = options.lowerbounds;
    tmp = sortrows([dose,response]);                                        % sort data in order of increasing dose
    opts.StartPoint = [options.baseline, options.upperbounds(2)./2, mean(options.xrange), max(response)];% min expected EMG noise; guess; 3/4 up the intensity scale; max MEP
    opts.Upper = options.upperbounds;                                       % max expected EMG noise; guess; twice max intensity; twice max M-wave
    [xData, yData] = prepareCurveData(tmp(:,1), tmp(:,2));                  % not sure what this does...

    % fit model to data____________________________________________________
    [fitresult, gof] = fit(xData, yData, ft, opts);                         % fit the data
    yfit = fitresult.d + (fitresult.a - fitresult.d) ./ (1 + (xData ./ fitresult.c).^fitresult.b );% compute the expected y-values

    % plot the fit_________________________________________________________
    if options.plot
        offset(1) = abs(max(dose) - min(dose)) .* 0.017;                    % shift text by 1.7% of the x range
        offset(2) = abs(max(response) - min(response)) .* 0.033;            % shift text by 3.3% of the y range
        plot([options.maxrange],[fitresult.a,fitresult.a],'r--');           % plot fitted baseline
        plot([options.maxrange],[fitresult.d,fitresult.d],'r--');           % plot fitted plateau
        plot([fitresult.c,fitresult.c],[fitresult.a,mean([fitresult.a,fitresult.d])],'r--');% plot fitted halfway point
        plot([options.maxrange(1),fitresult.c],[mean([fitresult.a,fitresult.d]),mean([fitresult.a,fitresult.d])],'r--');% plot to the fitted halfway point
        plot(tmp(:,1),yfit,'r-','LineWidth',1.5);                           % add fit to the plot
        a=axis;
        text(options.maxrange(1),a(4).*0.950, ['  y = ',num2str(fitresult.d,3),' + (',num2str(fitresult.a,3),'-',num2str(fitresult.d,3),') / (1+ (x / ',num2str(fitresult.c,3),')^{',num2str(fitresult.b,3),'}, r^2 = ',num2str(gof.adjrsquare,3)],'Color','r');
	text(options.maxrange(1),a(4).*0.850, ['  plateau = ',num2str(fitresult.d,3)], 'Color','k','FontSize', 8);% plateau point on Y-axis
	text(options.maxrange(1),a(4).*0.775, ['  baseline = ',num2str(fitresult.a,3)], 'Color','k','FontSize', 8);% baseline point on Y-axis
	text(options.maxrange(1),a(4).*0.700, ['  growth rate = ',num2str(fitresult.b,3)], 'Color','k','FontSize', 8);% baseline point on Y-axis
	text(options.maxrange(1),a(4).*0.625, ['  50% response = (',num2str(fitresult.c,3),', ',num2str(mean([fitresult.a,fitresult.d]),3),')'], 'Color','k','FontSize', 8);% 50% points on X and Y axes
    end
end