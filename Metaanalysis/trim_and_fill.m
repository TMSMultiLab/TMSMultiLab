%
% trim and fill algorithm for meta-analysis:
%
% [filled, stats] = trim_and_fill(x, s [, side] [, estimator] [, method] [, verbose])
%
% x = mx1 array of study means [required]
%
% s = mx1 array of study SEs [required]
%
% side = -1 for left, negative effects suppressed; 1 for right, positive effects suppressed (default = -1)
%
% estimator = 'L0', 'R0', 'Q0' (default = 'L0')
%
% method = 'FE','DL','HM','HO','HS','SJ' (default = 'DL' Dersimonian-Laird)
%
% verbose = false, true (default = false)
%
% from Taylor & Tweedle (2000)?
% There are n studies trying to measure effect D. Each study, j, produces an effect size Yj, which estimates D with variance Sj
% In addition to n observed studies, there are k0 missing
% suppose that there are originally n + k0 studies with this structure, and that k0 values from the set Xj have been suppressed, leaving a set of n observed studies.
% 'Our key assumption is: The suppression has taken place in such a way that it is the k0 values of the Xj with the most extreme negative ranks that have been suppressed
%
% The non-suppressed studies, n, are re-ranked. Three statistics are generated from this set of ranks:
% gamma *> 0 = the right-most run of positive ranks
% Tn = Wilcoxon statistic for the n ranked studies, then
% R0 = gamma*-1
% L0 = [4Tn - n(n+1)] / (2n+1)
% Q0 = n - 1/2 - sqrt(2n^2 - 4Tn + 1/4)
%
% Right = positive effects (most extreme effects)
%
% Algorithm:
% 1. Estimate meta-analytic mean effect, D, get initial centred values, and estimate of k0
% 2. Remove right-most (extreme) values, and re-estimate D using the trimmed set of data; construct centred values and new k0(2) from this set
% 3. Remove k0(2) values from right-most end, re-estimate D using trimmed data
% 4. Continue until the number of removed points, k0, is the same across 2 iterations
% 5. Fill the funnel plot with all the k0 symmetric points
%
% REFERENCES:
% https://pmc.ncbi.nlm.nih.gov/articles/PMC6571372/
% https://www.rdocumentation.org/packages/metafor/versions/2.4-0/topics/trimfill
% https://github.com/wviechtb/metafor/blob/master/R/trimfill.rma.uni.r

function [filled, stats] = trim_and_fill(x, s, side, estimator, method, verbose)

    %% PROCESS THE INPUT ARGUMENTS_________________________________________
    switch nargin
        case 0
            error('trimfill: no data specified');
        case 1
            error('trimfill: no standard error specified');
        case 2
            side = [];
            estimator = 'L0';
	    method = 'DL';
	    verbose = false;
        case 3
            estimator = 'L0'; 
            method = 'DL';	
            verbose = false;
        case 4
            method = 'DL';
	    verbose = false;
        case 5
	    verbose = false;
    end
    if size(x,2)>size(x,1)
        x = x';
    end
    if size(x,2)>2
        warning('trimfill: matrix entered as data, when m x 1 array needed: means');
    end
    if numel(x) < 3
        error('trimfill: too few datapoints');
    end
    if size(s,2)>size(s,1)
        s = s';
    end
    if size(s,2)>2
        warning('trimfill: matrix entered as data, when m x 1 array needed: SEs');
    end
    if numel(s) < 3
        error('trimfill: too few datapoints: SEs');
    end
    if numel(x) ~= numel(s)
        error('trimfill: different number of mean and SE datapoints entered');
    end
    if side ~= -1 & side ~= 1
        side = [];
    end
    if sum(strcmp(estimator, {'L0','RO','Q0'})) == 0
        estimator = 'L0';
    end
    if sum(strcmp(method, {'FE','DL','HM','HO','HS','SJ'})) == 0
        method = 'DL';
    end
    if ~verbose
        verbose = false;
    end

    %% FLIP THE DATA ?_____________________________________________________
    if side == 1                                                            % if side is POSITIVE
        x = -x;                                                             % flip the data
    elseif isempty(side)                                                    % if side not specified
        if max(x) < max(abs(x))                                             % use the least-extreme side
            side = 1;                                                       % set this to the positive side
            x = -x;                                                         % flip the data
	else
	    side = -1;                                                      % set this to the negative side
        end
    end


    %% SORT THE DATA_______________________________________________________
    [x, ix] = sortrows(x);                                                  % sort means by increasing size
    s = s(ix);                                                              % sort SEs by increasing means
    k = numel(x);                                                           % how many studies entered
    k0 = -1;                                                                % set the N suppressed studies to -1
    k0_tmp = 0;                                                             % current estimate of N suppressed studies
    iterations = 10;                                                        % how many times to try?
    i = 0;                                                                  % current iteration


    %% START THE ALGORITHM_________________________________________________
    while abs(k0_tmp - k0) > 0                                              % while the current estimate is different from the previous
        k0 = k0_tmp;                                                        % reset k0
        i = i+1;                                                            % increment iterations
        if i>iterations
            warning(['trimfill: algorithm did not converge in ',int2str(iterations),' iterations']);% warning message
	    stats.converge = false;
            break;                                                          % exit the while loop
        end

        % truncate the data by removing the latest estimet of k0 suppressed studies
        x_truncated = x(1:k-k0_tmp);                                        % remove the most-extreme k0 values off the end of means
        s_truncated = s(1:k-k0_tmp);                                        % remove the most-extreme k0 values off the end of SEs

        % get new estimate of meta effect size using truncated data________
	meta = meta_analysis(x_truncated, s_truncated, .05, 0, method);
	if method == 'FE'
            M = meta.FE_Mean;
        else
	    M = meta.RE_Mean;
        end
        % centre, rank and sign the original data around the new truncated meta effect size
        x_cent = x - M;                                                     % centre the values around the new meta-mean estimate
        [~,~,x_ranks] = unique(abs(x_cent));                                % rank the absolute centred values
        

        %% ADJUST RANKS FOR TIES___________________________________________
        for r=1:k                                                           % for each possible rank
            ix = x_ranks==r;                                                % which ranks are the same?
            if  sum(ix) > 1                                                 % if there are ties
                ix2 = x_ranks > r;                                          % which remaining ranks are bigger than the tied one?
                x_ranks(ix2) = x_ranks(ix2) + sum(ix)-1;                    % add ties-1 to each remaining rank
                x_ranks(ix) = mean(r:r+sum(ix)-1);                          % reset the tied ranks to the mean of ranks
                r = r+sum(ix)-1;                                            % move forward to the next rank
            end
        end
        x_signs = x_cent./abs(x_cent);                                      % signs of the ranked centred values


        %% ESTIMATE NUMBER OF SUPPRESSED STUDIES___________________________

        % R0 estimator_____________________________________________________
        gamma = k - max(-1.*(x_ranks.*x_signs).*((x_ranks.*x_signs)<0));    % gamma*>0 = the right-most run of positive ranks
        R0 = gamma - 1;                                                     % excess right-run positive ranks
        R0_se = sqrt(2.*max(0,R0) + 2);                                     % SE for this

        % Wilcoxon statistic for the n ranked studies______________________
        Tn = sum(x_ranks.*x_signs);                                         % sum of signed ranks = Tn

        % L0 estimator_____________________________________________________
	Sr = sum((x_ranks.*x_signs).*((x_ranks.*x_signs)>0));               % sum of positive-signed ranks
        L0 = ((4.*Sr) - (k.*(k + 1))) ./ (2.*k - 1);
        L0_Sr_var = (1./24) .* (k.*(k+1).*(2.*k+1) + 10.*L0.^3 + 27.*L0.^2 + 17.*L0 - 18.*k.*L0.^2 - 18.*k.*L0 + 6.*k.^2.*L0);
        L0_se = 4.*sqrt(L0_Sr_var) ./ (2.*k - 1);

        % Q0 estimator
        Q0 = k - 1/2 - sqrt(2.*k.^2 - 4.*Tn + 0.25);
        Q0_Tn_var = (1./24) .* (k.*(k+1).*(2.*k+1) + 10.*Q0.^3 + 27.*Q0.^2 + 17.*Q0 - 18.*k.*Q0.^2 - 18.*k.*Q0 + 6.*k.^2.*Q0);
        Q0_se = 2.*sqrt(Q0_Tn_var) ./ sqrt((k-0.5).^2 - k0_tmp.*(2.*k - k0_tmp - 1));


        %% SET UP NEXT ITERATION___________________________________________
        switch estimator
            case 'R0'
                k0_tmp = R0;
            case 'L0'
                k0_tmp = L0;
            case 'Q0'
                k0_tmp = Q0;
        end

        % round K0 & make it non-negative__________________________________
        k0_tmp = max(0, round(k0_tmp));
        % se.k0 <- max(0, se.k0)

        if verbose
	    disp([' iteration: ',int2str(i),', ',method,', mean = ',num2str(M,3),', k0 = ',int2str(k0_tmp),'...']);
	end
    end

    %% OUTPUT THE STATS____________________________________________________
    if i<iterations
        stats.converge = true;
    end
    stats.estimator = estimator;
    stats.side = side;
    stats.method = method;
    stats.i = i;
    stats.x_cent = x_cent;
    stats.x_ranks = x_ranks;
    stats.x_signs = x_signs;
    stats.k = k;
    stats.k0 = k0;
    stats.k0_tmp = k0_tmp;
    stats.L0 = L0;
    stats.R0 = R0;
    stats.Q0 = Q0;
    stats.Tn = Tn;
    stats.M = M;
    
    %% FILL IN THE SUPRESSED STUDIES ?_____________________________________
    if k0_tmp==0
        filled.x = x;
        filled.s = s;
    else
        filled.x = [-x_cent(k:-1:k-k0+1);x_cent];
        filled.s = [s(k:-1:k-k0+1);s];
	meta = meta_analysis(x_truncated, s_truncated, .05, 0, method);
	if method == 'FE'
            filled.M = meta.FE_Mean;
        else
	   filled. M = meta.RE_Mean;
        end
    end

    %% UNFLIP THE DATA_____________________________________________________
    if side ==1
        filled.x = -filled.x;                                               % flip the data back to the original direction
	stats.M = -stats.M;
    end
end