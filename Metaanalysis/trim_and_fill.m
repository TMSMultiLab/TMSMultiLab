%% trim and fill algorithm for meta-analysis

% from Taylor & Tweedle (2000)?
% There are n studies trying to measure effect D. Each study, j, produces an effect size Yj, which estimates D with variance Sj
% In addition to n observed studies, there are k0 missing
% suppose that there are originally n + k0 studies with this structure, and that k0 values from the set Xj have been suppressed, leaving a set of n observed studies.
% "Our key assumption is: The suppression has taken place in such a way that it is the k0 values of the Xj with the most extreme negative ranks that have been suppressed

% The non-suppressed studies, n, are re-ranked. Three statistics are generated from this set of ranks:
% gamma*>0 = the right-most run of positive ranks
% Tn = Wilcoxon statistic for the n ranked studies, then
% R0=gamma*-1
% L0=[4Tn - n(n+1)]/(2n+1)
% Q0=n - 1/2 - sqrt(2n^2 - 4Tn + 1/4)

% Right = positive effects (most extreme effects)
%1. Estimate D, get initial centred values, and estimate of k0
%2. Remove right-most values, and re-estimate D using the trimmed set of data; construct centred values and new k0(2) from this set
%3. Remove k0(2) values from right-most end, re-estimate D using trimmed data
%4. Continue until the number of removed points, k0, is the same across 2 iterations
%5. Fill the funnel plot with all the k0 symmetric points

function [filled, L0, R0, Q0, Tn] = trim_and_fill(x, s, side, estimator)

    %% PROCESS THE INPUT ARGUMENTS_________________________________________
    if nargin == 0
        error('trimfill: no data specified');
    elseif nargin == 1
        error('trimfill: no standard error specified');
    elseif nargin == 2
        side = [];
        estimator = "L0";
    elseif nargin ==3
        estimator = "L0";        
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


    %% START THE ANALYSIS__________________________________________________


    %% FLIP THE DATA ?_____________________________________________________
    if side == 1                                                            % if side is POSITIVE
        x = -x;                                                             % flip the data
    elseif isempty(side)                                                      % if side not specified
        if max(x) < max(abs(x))                                             % use the least-extreme side
            side = 1;                                                       % set this to the positive side
            x = -x;                                                         % flip the data
        end
    end


    %% SORT THE DATA_______________________________________________________
    [x, ix] = sortrows(x);                                                  % sort means by increasing size
    s = s(ix);                                                              % sort SEs by increasing means
    k = numel(x);                                                           % how many studies entered
    k0 = -1;                                                                % set the N suppressed studies to -1
    k0_tmp = 0;                                                             % current estimate of N suppressed studies
    iterations = 100;                                                       % how many times to try?
    i = 0;                                                                  % current iteration


    %% START THE ALGORITHM_________________________________________________
    while abs(k0_tmp - k0) > 0                                              % while the current estimate is different from the previous
        k0 = k0_tmp;                                                        % reset k0
        i = i+1;                                                            % increment iterations
        if i>iterations
            warning(['trimfill: algorithm did not converge in ',int2str(iterations),' iterations']);% warning message
            break;                                                          % exit the while loop
        end

        % truncate the data by removing the latest estimet of k0 suppressed studies
        x_truncated = x(1:k-k0);                                            % remove the most-extreme k0 values off the end of means
        s_truncated = s(1:k-k0);                                            % remove the most-extreme k0 values off the end of SEs

        % get new estimate of meta effect size using TRUNCATED data %%% DO THIS!!! %%%
        M = 0;

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
        L0 = ((4.*Tn) - (k.*(k + 1))) ./ (2.*k - 1);
        L0_Tn_var = (1./24) .* (k.*(k+1).*(2.*k+1) + 10.*L0.^3 + 27.*L0.^2 + 17.*L0 - 18.*k.*L0.^2 - 18.*k.*L0 + 6.*k.^2.*L0);
        L0_se = 4.*sqrt(L0_Tn_var) ./ (2.*k - 1);

        % Q0 estimator
        Q0 = k - 1/2 - sqrt(2.*k.^2 - 4.*Tn + 0.25);
        Q0_Tn_var = (1./24) .* (k.*(k+1).*(2.*k+1) + 10.*Q0.^3 + 27.*Q0.^2 + 17.*Q0 - 18.*k.*Q0.^2 - 18.*k.*Q0 + 6.*k.^2.*Q0);
        Q0_se = 2.*sqrt(Q0_Tn_var) ./ sqrt((k-0.5).^2 - k0.*(2.*k - k0 - 1));


        %% SET UP NEXT ITERATION___________________________________________
        switch estimator
            case "R0"
                k0_tmp = R0;
            case "L0"
                k0_tmp = L0;
            case "Q0"
                k0_tmp = Q0;
        end

        % round K0 & make it non-negative__________________________________
        k0 = max(0, round(k0));
        % se.k0 <- max(0, se.k0)

    end

    %% FILL IN THE SUPRESSED STUDIES ?_____________________________________
    if k0==0
        filled.x = x;
        filled.s = s;
    else
        filled.x = x;
        filled.s = s;
    end

    %% UNFLIP THE DATA_____________________________________________________
    if side ==1
        filled.x = -filled.x;                                               % flip the data back to the original direction
    end
end

% SEE:
% https://pmc.ncbi.nlm.nih.gov/articles/PMC6571372/
% https://www.rdocumentation.org/packages/metafor/versions/2.4-0/topics/trimfill
% https://github.com/wviechtb/metafor/blob/master/R/trimfill.rma.uni.r

%    while (abs(k0 - k0.sav) > 0) {
% 
%    #########################################################################
% 
%    ### if estimated number of missing studies is > 0
% 
%    if (k0 > 0) {
% 
%       ### flip data back if side is right
% 
%       if (side == "right") {
%          yi.c <- -1 * (yi.c - beta)
%       } else {
%          yi.c <- yi.c - beta
%       }
% 
%       ### create filled-in data set
% 
%       yi.fill <- c(x$yi.f, -1*yi.c[(k-k0+1):k])
% 
%       ### apply limits if specified
% 
%       if (!missing(ilim)) {
%          ilim <- sort(ilim)
%          if (length(ilim) != 2L)
%             stop(mstyle$stop("Argument 'ilim' must be of length 2."))
%          yi.fill[yi.fill < ilim[1]] <- ilim[1]
%          yi.fill[yi.fill > ilim[2]] <- ilim[2]
%       }
% 
%       vi.fill <- c(x$vi.f, vi[(k-k0+1):k])
%       wi.fill <- c(x$weights.f, wi[(k-k0+1):k])
%       ni.fill <- c(x$ni.f, ni[(k-k0+1):k])
% 
%       ### add measure attribute to the yi.fill vector
% 
%       attr(yi.fill, "measure") <- x$measure
% 
%       ### fit model with imputed data
% 
%       args <- list(yi=yi.fill, vi=vi.fill, weights=wi.fill, ni=ni.fill, method=x$method, weighted=x$weighted, digits=x$digits, ...)
%       res <- suppressWarnings(.do.call(rma.uni, args))
% 
%       ### fill, ids, and slab are of length 'k.f + k0' (i.e., subsetted but with NAs)
% 
%       res$fill <- c(rep(FALSE,x$k.f), rep(TRUE,k0))
%       res$ids  <- c(x$ids, (max(x$ids)+1):(max(x$ids)+k0))
% 
%       if (x$slab.null) {
%          res$slab <- c(paste("Study", x$ids), paste("Filled", seq_len(k0)))
%       } else {
%          res$slab <- c(x$slab, paste("Filled", seq_len(k0)))
%       }
%       res$slab.null <- FALSE
% 
%    } else {
% 
%       ### in case 0 studies are imputed
% 
%       res      <- x
%       res$fill <- rep(FALSE,k)
% 
%    }
% 
%    res$k0     <- k0
%    res$se.k0  <- se.k0
%    res$side   <- side
%    res$k0.est <- estimator
%    res$k.all  <- x$k.all + k0
% 
%    if (estimator == "R0") {
%       m <- -1:(k0-1)
%       res$p.k0 <- 1 - sum(choose(0+m+1, m+1) * 0.5^(0+m+2))
%    } else {
%       res$p.k0 <- NA_real_
%    }
% 
%    class(res) <- c("rma.uni.trimfill", class(res))
%    return(res)
% }