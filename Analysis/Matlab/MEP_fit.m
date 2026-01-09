%% [fitresult, gof, options] = MEP_fit(timepoints, data [,options])
% fit a non-linear model to an MEP
%  INPUTS
%     timepoints: m x 1 array of timepoints in ms, eg 10:0.25:50)
%     data: m x 1 array of EMG data, in mV
%  Output:
%      fitresult : a fit object representing the fit
%      gof : structure with goodness-of fit info
%
% FITTING SHOULD OPTIMISE:
%     latency (muscle-dependent, typical range of ~18-25 for hand muscles; min/max = 5,60?)
%     offset (DC offset, typically 0mV, min/max +- 1mV)
%     scale (MEP amplitude, typically ~0.2 to give 1mV; min/max = +- 5 to give ~20mV MEP)
%     asymmetry (typical range of ~ +-3; min/max = +- 10)
%
% the full equation for each component is:
%     scale1 .* (-((x-offset_x1) + asymmetry1).*exp(-((x-offset_x1)/lambda1).^2)) + offset_y1 + 
%     scale2 .* (-((x-offset_x2) + asymmetry2).*exp(-((x-offset_x2)/lambda2).^2)) + offset_y2 +
%     ...                                                                                     +
%     scalen .* (-((x-offset_xn) + asymmetryn).*exp(-((x-offset_xn)/lambdan).^2)) + offset_yn;

function [fitresult, gof, output] = MEP_fit(timepoints, data, options)

    %% OPTIONS ... waiting for future development
    %options.plot = true;
    
    output = [];

    %% Fit: MEP fit________________________________________________________
    [xData, yData] = prepareCurveData( timepoints, data );

    %% SET UP FIT TYPE AND OPTIONS_________________________________________
    % fit 2x MUAP models to the MEP
    % a: scale = multiplier to fit the peak-to-peak amplitude
    % b: offset_x = subtracted from the timepoints to scale the model around 0 => latency
    % c: asymmetry = how different are the two phases of the biphasic MUAP shape? -10:10, 0 = symmetrical
    % d: lambda = width of the function
    % e: offset_y = DC offset on the y-axis - ideally set to ~zero by de-meaning the pre-MEP data

    ft = fittype( 'a1.*(-((x-b1)+c1).*exp(-((x-b1)/d1).^2))+e1 + a2.*(-((x-b2)+c2).*exp(-((x-b2)/d2).^2))+e2', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    %                   a1     a2    b1 b2 c1 c2  d1 d2 e1 e2
    opts.Lower =      [-10    -10    20 20 -1 -1  0   0 -0.25 -0.25];       % lower bounds
    opts.StartPoint = [  0.25   0.25 30 30  0  0  3   3  0     0];          % starting point
    opts.Upper =      [ 10     10    40 40  1  1 10  10  0.25  0.25];       % upper bounds
    % FREE               *      *     *  *  *  *  *   *  *     *

    % Fit model to data.
    [fitresult, gof] = fit(xData,yData,ft,opts);
    
    y1 = fitresult.a1 .* (-((timepoints-fitresult.b1) + fitresult.c1).*exp(-((timepoints-fitresult.b1)/fitresult.d1).^2)) + fitresult.e1;% compute the fit
    output.auc1 = sum(abs(y1-fitresult.e1));                                % area under the curve (after subtracting y-offset)

    y2 = fitresult.a2 .* (-((timepoints-fitresult.b2) + fitresult.c2).*exp(-((timepoints-fitresult.b2)/fitresult.d2).^2)) + fitresult.e2;% compute the fit
    output.auc2 = sum(abs(y2-fitresult.e2));                                % area under the curve (after subtracting y-offset)

    % compute total AUC under fit and model________________________________
    output.auc_fit = output.auc1 + output.auc2;                             % area under the modelled components may be > auc the data
    output.auc_data = sum(abs(data));
    
    if options.plot
        % Create a figure for the plots
        figure( 'Name', 'MEP fit' );

        %% PLOT DATA + FIT_____________________________________________________
        subplot(3,1,1);
        hold on;
        grid on;
        plot( fitresult, xData, yData );
        legend off;
        ylabel( 'EMG, mV', 'Interpreter', 'none' );                             % Label axes
        text(min(xData)+(max(xData)-min(xData)).*0.025,max(yData),['R^2 = ',num2str(gof.adjrsquare,3)]);% model fit
        title('Data + fit');
        a=axis;

        %% PLOT RESIDUALS______________________________________________________
        subplot(3,1,2);
        hold on;
        grid on;
        plot( fitresult, xData, yData, 'residuals' );
        legend off;
        ylabel( 'EMG, mV', 'Interpreter', 'none' );                             % Label axes
        title('Residuals');

        %% PLOT INDIVIDUAL FIT COMPONENTS______________________________________
        subplot(3,1,3);
        hold on;
        grid on;
        plot(timepoints,y1,'r-');                                               % plot component 1 in red
        plot(timepoints,y2,'b-');                                               % plot component 1 in blue
	
	%% FORMAT THE PLOT_____________________________________________________
        xlabel( 'Time after stimulus, ms', 'Interpreter', 'none' );             % Label axes
        ylabel( 'EMG, mV', 'Interpreter', 'none' );                             % Label axes
        text(min(xData)+(max(xData)-min(xData)).*0.025,max(yData),['AUC model = ',num2str(output.auc_fit,3),'; data = ',num2str(output.auc_data,3)]);% AUC
        title('Components');
        axis(a);
        set(gcf,'Position',[0,40,600,600]);
    end
end