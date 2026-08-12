%% CREATE A FOREST PLOT FROM A LIST OF MEANS, SEs (, study labels and options)
% function [options, stats] = forest_plot(Ms, SEs[, labels] [, options])
% Ms = mx1 array of study means
% SEs = mx1 array of study standard errors
% labels = mx1 array of study labels as strings
% opts = structure of options, any of which may be specified or not, including:
%   opts.subgroups = 0; default = no subgroups. If one integer per mean provided, these are used to sub-group the data
%   opts.order = default = order of input (1 value per mean; unique integers)
%   opts.sort = 'default' ['ascending', 'descending', 'custom'] - plot the data in this order; default = order of input, custom requires opts.order
%   opts.null = 0 - what counts as zero effect?
%   opts.alpha = 0.05 - alpha level for studywise confidence intervals
%   opts.symbol = 's' - symbol for study means (can be mx1 array of strings (not cells), eg: ['s';'o';'d';...])
%   opts.linestyle = '-' - linestyle for study CIs (can be mx2 array of strings (not cells), eg: ['- ';': ';'--';...])
%   opts.markersize = 12 - how large the mean symbol size (adjusted by meta-analysis weights to 0.4:2* this)
%   opts.color = 'k' - colour for study lines and symbols (can be RGB triplet, mx1 array of strings (not cells) or mx3 vector of triplets)
%       eg: 'k';  [0,1,0];  ['k';'r';'g';...];  [0,0,1;0,1,0;1,0,0;...]
%   opts.linewidth = 2 - linewidth for study lines
%   opts.fontsize = [16,16] size of fonts for Figure and Axes (modified for large N)
%   opts.axisposition = [0.2, 0.1, 0.55, 0.8]- default size of the main plot, allowing for text labels
%   opts.figureposition = [1, 50, 800, 800] - default figure size
%
% outputs opts (as above, some modified, some new) and stats - the meta-analysis stats (see meta_analysis function for details)

function [opts, stats] = forest_plot(Ms, SEs, labels, opts)

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
        opts.color = 'k';                                                   % initialise options here, then complete below...
    end

    %% IF SOME OF THE OPTIONS ARE NOT SPECIFIED____________________________
    if ~isfield(opts,'subgroups')
        opts.subgroups = ones(numel(Ms),1);                                 % default subgroups = mx1
        opts.subgrouplabels = '';                                           % empty subgroup labels
    end
    if ~isfield(opts,'subgrouplabels')
        opts.subgrouplabels = '';                                           % empty subgroup labels
    end
    if numel(opts.subgroups) ~= numel(Ms)
        warning('forest_plot: wrong number of subgroups enteredl, defaulting to none');
        opts.subgroups = ones(numel(Ms),1);                                 % default subgroups = mx1
    end
    if ~isfield(opts,'sort')
        opts.sort = 'default';                                              % default sort = 1:k
    end 
    if ~isfield(opts,'null')
        opts.null = 0;                                                      % default null = 0
    end
    if ~isfield(opts,'alpha')
        opts.alpha = 0.05;                                                  % default alpha = 0.05
    end
    if ~isfield(opts,'symbol')
        opts.symbol = 's';                                                  % symbol for studies
    end
    if ~isfield(opts,'markersize')
        opts.markersize = 12;                                               % default (mean) MarkerSize property
    end
    if ~isfield(opts,'linestyle')
        opts.linestyle = '-';                                               % linestyle for studies
    end
    if ~isfield(opts,'color')
        opts.color = 'k';                                                   % linecolour for studies
    end
    if ~isfield(opts,'linewidth')
        opts.linewidth = 2;                                                 % linewidth for studies
    end
    if ~isfield(opts,'fontsize')
        opts.fontsize = [16, 16];                                           % fontsize for [Figure, Axes] (Labels = 0.85*)
    end   
    if ~isfield(opts,'axisposition')
        opts.axisposition = [0.24, 0.1, 0.54, 0.8];                         % size of the plot (allows for text on left & right)
    end 
    if ~isfield(opts,'figureposition')
        opts.figureposition = [1, 51, 750, 800];                            % default figure size
    end    
    if ~isfield(opts,'studyoffset')
        opts.studyoffset = 0.42;                                            % default study text offset proportion
    end  
    if ~isfield(opts,'sfs')
        opts.sfs = 3;                                                       % significant figures in the text labels
    end
    if ~isfield(opts,'scalek')
        opts.scalek = 13;                                                   % if more than this number of studies, decrease fonts
    end
    subgroups = unique(opts.subgroups);                                     % ordered list of subgroups
    if size(opts.subgroups,2)>size(opts.subgroups,1)
        opts.subgroups = opts.subgroups';                                   % rotate the subgroups
    else
        subgroups = subgroups';                                             % rotate the index to the subgroups
    end
    if numel(opts.subgrouplabels) ~= subgroups
        warning('forest_plot: wrong number of subgroup labels specified, defaulting to none');
        opts.subgrouplabels = '';                                           % empty subgroup labels
    end
    
    %% ORDER THE DATA______________________________________________________
    if numel(subgroups) == 1                                                % if only one subgroup specified
        switch opts.sort
            case 'default'
                jx = 1:numel(Ms);                                           % no sorting - order of arrival
        case 'ascending'
            [~, jx] = sortrows(Ms,'ascend');                            % ascending order of means
        case 'descending'
            [~, jx] = sortrows(Ms,'decend');                            % descending order of means
        case 'custom'
            if ~isfield(opts,'order')
                warning('forest_plot: no sort order specified, defaulting to input');
                jx = 1:numel(Ms);
            else
                if numel(opts.order) ~= numel(Ms)
                    warning('forest_plot: sort order different length to data, defaulting to input');
                    jx = 1:numel(Ms);
                else
                    jx = opts.order;
                end
            end
        end
    else
        [~, jx] = sortrows(opts.subgroups,'ascend');                        % ascending order of subgroups (no options here yet)
    end
    Ms = Ms(jx);
    SEs = SEs(jx);
    labels = labels(jx);
    if size(opts.color,1)==numel(Ms)
        opts.color = opts.color(jx,:);
    end


    %% INTERIM NUMBER OF STUDIES___________________________________________
    k = numel(Ms);                                                          % number of studies
    if k>opts.scalek                                                        % for k>opts.scalek, adjust plot sizes...
        opts.markersize = opts.markersize  .* opts.scalek./k;               % adjust mean size of markers
        opts.linewidth = opts.linewidth    .* opts.scalek./k;               % adjust linewidth size
        opts.fontsize(2) = opts.fontsize(2).* opts.scalek./k;               % adjust fontsize for studies
    end


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
        opts.color = opts.color([ix;true],:);                               % new list of colours + the last one for the mean
        opts.symbol = opts.color(ix,:);                                     % new list of symbols
        opts.linestyle = opts.linestyle(ix,:);                              % new list of linestyles
        if numel(opts.subgroups)==numel(Mx);
            opts.subgroups = opts.subgroups(ix);                            % new list of subgroups
        end
    end


    %% COMPUTE CONSTANTS___________________________________________________
    k = numel(Ms);                                                          % number of studies
    if numel(subgroups)==1                                                  % extra lines to add for gaps between subgroups & after last one
        l = 0.5;                                                            % 1 group = 0.5
    else
        l = numel(subgroups)+0.5;                                           % 2 groups = 2.5, 3 groups = 3.5, etc
    end
    opts.conf = 100.*(1-opts.alpha);                                        % calculate confidence level


    %% COMPUTE CONFIDENCE INTERVALS________________________________________
    CIs = SEs .* norminv(1-(opts.alpha./2));
    opts.xlims = [min(Ms) - max(CIs).*1.05, max(Ms) + max(CIs).*1.05];      % x-axis limits


    %% META-ANALYSE THE FULL DATASET_______________________________________
    stats = meta_analysis(Ms, SEs, opts.alpha, 0, 'DL');                    % meta-analysis using Dersimonian-Laired method
    M_mean = stats.RE_Mean;                                                 % meta-analysis mean to plot
    M_CI = stats.RE_SE .* norminv(1-(opts.alpha./2));                       % meta-analysis CI to plot
    opts.scale = 0.5 + (stats.Ws - min(stats.Ws)) ./ (max(stats.Ws)-min(stats.Ws)); % scale the SEs for marker size adjustments between 0.5 and 1.5 x the default
    fmt = ['%.',int2str(opts.sfs - ceil(log10(M_mean))),'f'];               % format for decimal places in number outputs?
    % 100-1000    3   0dp
    % 10-100      2   1dp
    % 1-10        1   2dp
    % 0.1-1       0   3dp
    % 0.01-0.1    -1  4dp


    %% OPEN A NEW FIGURE___________________________________________________
    figure();
    hold on;
    set(gcf, 'Color', 'w');
    set(gca,'FontSize', opts.fontsize(1));
    box on;
    plot([opts.null, opts.null], [-1, k + l + 1], 'k-');                    % plot zero effect
    set(gca,'Position', opts.axisposition);


    %% PLOT THE DATA_______________________________________________________
    for n = 1:k
        ypos = k + l - opts.subgroups(n) + 2 - n;                           % where on the y-axis is this datapoint?
    
        % CONFIDENCE INTERVALS
        plot([Ms(n)-CIs(n), Ms(n)+CIs(n)], [ypos, ypos], [opts.linestyle(n,:)], 'Color',opts.color(n,:), 'LineWidth', opts.linewidth);

        % MEANS
        plot(Ms(n), ypos, opts.symbol(n,:), 'Color', opts.color(n,:), 'MarkerSize', opts.markersize.*opts.scale(n), 'MarkerFaceColor', opts.color(n,:));

        % LABELS ON THE LEFT
        text(opts.xlims(1)-abs(diff(opts.xlims)).*opts.studyoffset, ypos, labels{n}, 'Color',opts.color(n,:),'FontSize',opts.fontsize(2).*0.85);

        % EFFECTS ON THE RIGHT
        text(opts.xlims(2)+abs(diff(opts.xlims)).*0.02, ypos, [sprintf(fmt,Ms(n)),' [',sprintf(fmt,Ms(n)-CIs(n)),', ',sprintf(fmt,Ms(n)+CIs(n)),']'], 'Color',opts.color(n,:), 'FontSize', opts.fontsize(2).*0.85);
    end
    
    % LEFT AND RIGHT COLUMN LABELS_________________________________________
    text(opts.xlims(1)-abs(diff(opts.xlims)).*opts.studyoffset, (k+l).*1.075, 'Study', 'FontAngle', 'Italic', 'FontSize', opts.fontsize(1).*0.85);
    text(opts.xlims(2)+abs(diff(opts.xlims)).*0.02, (k+l).*1.075, ['Mean ',char(177),' ',int2str(opts.conf),'% CI'], 'FontAngle', 'Italic', 'FontSize', opts.fontsize(1).*0.85);

    % SUBGROUP MEANS_______________________________________________________
    if numel(subgroups) > 1
        for s = subgroups                                                   % for each unique subgroup number
        idx = find(opts.subgroups==s);                                  % index to means in this subgroup
            ms = Ms(idx);                                                   % subgroup means
        ses = SEs(idx);                                                 % subgroup SEs
        tmp = meta_analysis(ms, ses, opts.alpha, 0, 'DL');              % meta-analysis using Dersimonian-Laired method
            m_mean = tmp.RE_Mean;                                           % meta-analysis mean to plot                
            m_CI = tmp.RE_SE .* norminv(1-(opts.alpha./2));                 % meta-analysis CI to plot
        top = k + l - s +2 - idx(1) + 0.5;                              % the top of this subgroup, for the shaded bar & line
        bottom = k + l - s +2 - idx(end) - 0.75;                        % the bottom of this subgroup, for the shaded bar & line
            patch([m_mean-m_CI, m_mean-m_CI, m_mean+m_CI, m_mean+m_CI], [bottom, top, top, bottom], opts.color(idx(1),:), 'FaceAlpha', 0.075, 'EdgeAlpha', 0.075);           % shaded bar
            plot([m_mean, m_mean], [bottom, top], '--', 'Color', opts.color(idx(1),:));                                                                                      % broken line
            patch([m_mean-m_CI, m_mean, m_mean+m_CI, m_mean], [bottom, bottom-(k.*0.01), bottom, bottom+(k.*0.01)], opts.color(idx(1),:), 'EdgeColor', opts.color(idx(1),:));% solid diamond
        if numel(opts.subgrouplabels)>1
            tmp = [opts.subgrouplabels{s},' mean'];
        else
            tmp = 'Mean';
        end
            text(opts.xlims(1)-abs(diff(opts.xlims)).*opts.studyoffset, bottom, tmp, 'FontWeight','Bold', 'Color', opts.color(idx(1),:),'FontSize', opts.fontsize(2).*0.85); % Mean label on the left, and data:
        text(opts.xlims(2)+abs(diff(opts.xlims)).*0.02, bottom, [sprintf(fmt,m_mean),' [',sprintf(fmt,m_mean-m_CI),', ',sprintf(fmt,m_mean+m_CI),']'], 'FontWeight', 'Bold', 'Color', opts.color(idx(1),:), 'FontSize', opts.fontsize(2).*0.85);
        end
    end

    % GRAND MEAN, CI, & LABELS_____________________________________________
    patch([M_mean-M_CI, M_mean-M_CI, M_mean+M_CI, M_mean+M_CI], [-k.*0.075, (k+l).*1.05, (k+l).*1.05, -k.*0.075], opts.color(n+1,:), 'FaceAlpha',0.075, 'EdgeAlpha', 0.075);% shaded bar
    plot([M_mean, M_mean], [-k.*0.075, (k+l).*1.05], '--', 'Color', opts.color(n+1,:));                                                                                     % dotted line
    patch([M_mean-M_CI, M_mean, M_mean+M_CI, M_mean], [0, -k.*0.02, 0, k.*0.02], opts.color(n+1,:), 'EdgeColor', opts.color(n+1,:));                                        % solid diamond
    if numel(opts.subgrouplabels)>1
        tmp = 'Overall mean';
    else
        tmp = 'Mean';
    end
    text(opts.xlims(1)-abs(diff(opts.xlims)).*opts.studyoffset, 0, tmp, 'FontWeight','Bold', 'Color', opts.color(n+1,:),'FontSize', opts.fontsize(1).*0.85);                % Mean label on the left, and data:
    text(opts.xlims(2)+abs(diff(opts.xlims)).*0.02, 0, [sprintf(fmt,M_mean),' [',sprintf(fmt,M_mean-M_CI),', ',sprintf(fmt,M_mean+M_CI),']'], 'FontWeight', 'Bold', 'Color', opts.color(n+1,:), 'FontSize', opts.fontsize(1).*0.85);


    %% FORMAT THE AXES_____________________________________________________
    axis([opts.xlims, -k.*0.075, (k+l).*1.05]);                             % rescale the axes
    xlabel('Effect size');                                                  % effect-size on X-axis
    yticks([]);                                                             % remove y-axis labels
    set(gcf, 'Position', opts.figureposition);                              % set the figure position

end