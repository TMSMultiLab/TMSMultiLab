%% CREATE A FOREST PLOT FROM A LIST OF MEANS, SEs (, study labels and options)
% function options = forest_plot(Ms, SEs[, labels] [, options])
% Ms = mx1 array of study means
% SEs = mx1 array of study standard errors
% labels = mx1 array of study labels as strings
% opts = structure of options, any of which may be specified or not, including:
%   opts.null = 0 - what counts as zero effect?
%   opts.alpha = 0.05 - alpha level for studywise confidence intervals
%   opts.symbol = 's' - symbol for study means (can be mx1 array of strings (not cells), eg: ['s';'o';'d';...])
%   opts.linestyle = '-' - linestyle for study CIs (can be mx2 array of strings (not cells), eg: ['- ';': ';'--';...])
%   opts.color = 'k' - colour for study lines and symbols (can be RGB triplet, mx1 array of strings (not cells) or mx3 vector of triplets)
%       eg: 'k';  [0,1,0];  ['k';'r';'g';...];  [0,0,1;0,1,0;1,0,0;...]
%   opts.linewidth = 2 - linewidth for study lines
%   opts.xlims = [minx, maxx] - x-axis limits

function opts = forest_plot(Ms, SEs, labels, opts)


    %% PROCESS THE INPUT ARGUMENTS_________________________________________
    if nargin==1
        error('forest_plot: Minimum of two mx1 non-empty arrays required'); % no SEs provided
    end
    if isempty(Ms) || isempty(SEs)
        error('forest_plot: Minimum of two mx1 non-empty arrays required'); % either Ms or SEs are empty
    end
    if ndims(Ms)>2 || ndims(SEs)>2
        error('forest_plot: Two mx1 non-empty arrays required');            % Ms or SEs are more than 2-dimensional
    end
    if size(Ms,1) < size(Ms,2)
        Ms = Ms';                                                           % rotate Ms to the longest dimension
    end
    if size(SEs,1) < size(SEs,2)
        SEs = SEs';                                                         % rotate SEs to the longest dimension
    end
    if numel(Ms) ~= numel(SEs)
        error('forest_plot: Two mx1 non-empty arrays required');            % different number of datapoints in Ms and SEs
    end
    if nargin<3 || isempty(labels)
        labels = append('Study ',string(1:numel(Ms)));                      % use numbers to ID the studies
    end
    if nargin<4 || isempty(opts)
        opts.null = 0;                                                      % default null = 0
        opts.alpha = 0.05;                                                  % default alpha = 0.05
        opts.symbol = 's';                                                  % symbol for studies
        opts.markersize = 12;                                               % default (mean) MarkerSize property
        opts.linestyle = '-';                                               % linestyle for studies
        opts.color = 'k';                                                   % linecolour for studies
        opts.linewidth = 2;                                                 % linewidth for studies
    end

    %% IF SOME OF THE OPTIONS ARE SPECIFIED________________________________


    %% FOR TESTING_________________________________________________________
    %opts.color = rand(numel(Ms),3);
    %opts.symbol = ['s';'o';'d';'p';'h';'^';'v';'<';'>'];
    %opts.linestyle = ['- ';': ';'-.';'--';'- ';': ';'-.';'--';'- '];


    %% INTERIM NUMBER OF STUDIES___________________________________________
    k = numel(Ms);                                                          % number of studies


    %% CREATE COLOUR VECTOR________________________________________________
    switch size(opts.color,1)
        case 1
            opts.color = repmat(opts.color,k+1,1);                          % extend colour to k+1 matrix
        case k
            if size(opts.color,2)>1
                opts.color = [opts.color;[0,0,0]];
            else
                opts.color = [opts.color;'k'];                              % add black as default for the mean
            end
        case k+1
        otherwise
            warning('forest_plot: wrong number of colors, defaulting to black');% remove non-finite datapoints
            opts.color = repmat('k',k+1,1);
    end


    %% CREATE SYMBOL VECTOR________________________________________________
    switch size(opts.symbol,1)
        case 1
            opts.symbol = repmat(opts.symbol,k,1);                          % extend symbol to k matrix
        case k
        otherwise
            warning('forest_plot: wrong number of symbols, defaulting to squares');% remove non-finite datapoints
            opts.symbol = repmat('s',k,1);
    end


    %% CREATE LINESTYLE VECTOR_____________________________________________
    switch size(opts.linestyle,1)
        case 1
            opts.linestyle = repmat(opts.linestyle,k,1);                    % extend symbol to k matrix
        case k
        otherwise
            warning('forest_plot: wrong number of linestyle, defaulting to solid');% remove non-finite datapoints
            opts.linestyle = repmat('- ',k,1);
    end


    %% HOW MANY FINITE STUDIES ARE THERE___________________________________
    ix = isfinite(Ms) & isfinite(SEs);                                      % which datapoints have finite M and SE?
    if sum(ix) < numel(Ms)
        warning('forest_plot: non-finite datapoints have been removed');    % remove non-finite datapoints
        Ms = Ms(ix);
        SEs = SEs(ix);
        labels = labels(ix);
        opts.color = opts.color([ix;true],:);                               % last one for the mean
    end


    %% COMPUTE CONSTANTS___________________________________________________
    k = numel(Ms);                                                          % number of studies
    opts.conf = 100.*(1-opts.alpha);                                        % calculate confidence level
    opts.scale = 1.5 - (SEs - min(SEs)) ./ (max(SEs)-min(SEs));             % scale the SEs for marker size adjustments


    %% COMPUTE CONFIDENCE INTERVALS________________________________________
    CIs = SEs .* norminv(1-(opts.alpha./2));
    opts.xlims = [min(Ms) - max(CIs).*1.05, max(Ms) + max(CIs).*1.05];      % x-axis limits


    %% META-ANALYSE THE DATA_______________________________________________
    stats = meta_analysis(Ms, SEs, opts.alpha, 0, 'DL');                    % meta-analysis using Dersimonian-Laired method
    M_mean = stats.RE_Mean;
    M_SE = stats.RE_SE;
    M_CI = M_SE .* norminv(1-(opts.alpha./2));

    fmt = ['%.',int2str(3 - ceil(log10(M_mean))),'f'];                      % format for decimal places in number outputs?
    % 100-1000    3   0dp
    % 10-100      2   1dp
    % 1-10        1   2dp
    % 0.1-1       0   3dp
    % 0.01-0.1    -1  4dp


    %% OPEN A NEW FIGURE___________________________________________________
    figure();
    hold on;
    set(gcf, 'Color', 'w');
    box on;
    plot([opts.null,opts.null],[-1,k+1],'k-');                              % plot zero effect
    set(gca,'Position',[0.2,0.1,0.55,0.8]);


    %% PLOT THE DATA_______________________________________________________
    for n = 1:k
        % CONFIDENCE INTERVALS
        plot([Ms(n)-CIs(n),Ms(n)+CIs(n)],[k+1-n,k+1-n],[opts.linestyle(n,:)],'Color',opts.color(n,:),'LineWidth',opts.linewidth);

        % MEANS
        plot(Ms(n),k+1-n,opts.symbol(n,:),'Color',opts.color(n,:),'MarkerSize',opts.markersize.*opts.scale(n),'MarkerFaceColor',opts.color(n,:));

        % LABELS ON THE LEFT
        text(opts.xlims(1)-abs(diff(opts.xlims)).*0.32, k+1-n, labels{n}, 'Color',opts.color(n,:));

        % EFFECTS ON THE RIGHT
        text(opts.xlims(2)+abs(diff(opts.xlims)).*0.02, k+1-n, [sprintf(fmt,Ms(n)),' [',sprintf(fmt,Ms(n)-CIs(n)),', ',sprintf(fmt,Ms(n)+CIs(n)),']'], 'Color',opts.color(n,:));
    end
    text(opts.xlims(1)-abs(diff(opts.xlims)).*0.32, k+1, 'Study','FontAngle','Italic');
    text(opts.xlims(2)+abs(diff(opts.xlims)).*0.02, k+1, ['Mean ',char(177),' ',int2str(opts.conf),'% CI'],'FontAngle','Italic');

    % MEANS & CI
    patch([M_mean-M_CI, M_mean, M_mean+M_CI, M_mean],[0, -0.25, 0, 0.25], opts.color(n+1,:), 'EdgeColor',opts.color(n+1,:));

    % LABEL ON THE LEFT 
    text(opts.xlims(1)-abs(diff(opts.xlims)).*0.32, 0, 'Mean', 'FontWeight','Bold', 'Color',opts.color(n+1,:));

    % EFFECT ON THE RIGHT
    text(opts.xlims(2)+abs(diff(opts.xlims)).*0.02, 0, [sprintf(fmt,M_mean),' [',sprintf(fmt,M_mean-M_CI),', ',sprintf(fmt,M_mean+M_CI),']'], 'FontWeight','Bold', 'Color',opts.color(n+1,:));


    %% FORMAT THE AXES_____________________________________________________
    axis([opts.xlims,-0.5,k+0.5]);                                          % rescale the axes
    xlabel('Effect size');                                                  % effect-size on X-axis
    yticks([]);                                                             % remove y-axis labels

end