function stats = meta_analysis(Ms, Vs, alpha, comparison, method)
% https://www.cochrane.org/authors/handbooks-and-manuals/handbook/current/chapter-10
% https://methods.cochrane.org/sites/methods.cochrane.org.statistics/files/uploads/SMG_training_course_2016/cochrane_smg_training_2016_viechtbauer.pdf
% https://jasp-stats.org/2017/11/15/meta-analysis-jasp/
% https://www.metafor-project.org/doku.php/metafor
% https://cjvanlissa.github.io/Doing-Meta-Analysis-in-R/heterogeneity-statistics.html
% https://www.jstatsoft.org/index.php/jss/article/view/v036i03/409
    % 1 = which estimator for Tau-squared?
    % 2 = weight the studies by their inverse variance

    %% PROCESS THE INPUTS________________________________________________
    switch nargin
        case 0
	    error('meta_analysis: requires at least two mx1 arrays: means and variances');
	    
	    case 1
	        error('meta_analysis: requires at least two mx1 arrays: means and variances');
	    
        case 2
            alpha = 0.05;
	        comparison = 0;
	        method = 'DL'; % default method
	
        case 3
            comparison = 0;
	        method = 'DL'; % default method
	
        case 4
    	    method = 'DL'; % default method
    end
    switch method
        case {'Fixed effects','fixed effects','fe','Fixed','fixed'}
	        method = 'FE';
        case {'DerSimonian-Laird','DSL','dsl','DL','dl','dersimonian-laird'}
	        method = 'DL';
	    case {'Positive DerSimonian and Laird','DLp'}
	        method = 'DL';
	        warning('meta_analysis: DLp method not implemented, using DL');
	    case {'Two-step DerSimonian and Laird','DL2'}
	        method = 'DL';
	        warning('meta_analysis: DL2 method not implemented, using DL');
	    case {'Hedges-Olkin','hedges-olkin','ho','Cochran','cochran'}
	        method = 'HO';
	    case {'Two-step Hedges and Olkin','HO2','ho2'}
	        method = 'HO';
	        warning('meta_analysis: HO2 method not implemented, using HO');
	    case {'Hartung-Makambi','HM'}
	        method = 'HM';
	        warning('meta_analysis: HM method not checked');
	    case {'Hunter-Schmidt','hs'}
	        method = 'HS';
	    case {'Paule-Mandel','PM','pm','generalized Q-statistic'}
	        method = 'PM';
	        error('meta_analysis: PM method not implemented');
	    case {'Maximum likelihood','maximum likelihood','ML','MLE','ml','mle'}
	        method = 'ML';
	        error('meta_analysis: ML method not implemented...');
	    case {'Restricted maximum likelihood','restricted maximum likelihood','REMLE','REML','ReML','ReMLE','reml','remle'}
	        method = 'REML';
	        error('meta_analysis: REML method not implemented...');
	    case {'Approximate restricted maximum likelihood', 'AREML','areml'}
	        method = 'AREML';
	        error('meta_analysis: AREML method not implemented...');
	    case {'Sidik-Jonkman','SJ','sj'}
	        method = 'SJ';
	    case {'Rukhin Bayes','RB'}
	        method = 'RB';
	        error('meta_analysis: RB method not implemented...');
	    case {'Positive Rukhin Bayes','RBp'}
	        method = 'RBp';
	        error('meta_analysis: RBp method not implemented...');
	    case {'Full Bayes'}
	        method = 'FB';
	        error('meta_analysis: FB method not implemented...');
	    case {'Bayes modal'}
	        method = 'BM';
	        error('meta_analysis: BM method not implemented...');
	    case {'Empirical Bayes'}
            method = 'EB';
	        error('meta_analysis: EB method not implemented...');
        otherwise
	        error('meta_analysis: method not recognised...');
    end
    if isempty(alpha) || alpha < 0 || alpha > 1
        alpha = 0.05;
    end
    if isempty(comparison | ~isnumeric(comparison))
        comparison = 0;
    end


    %% COMMON STATS NEEDED_______________________________________________
    stats.method = method;
    stats.Alpha = alpha;
    stats.null = comparison;
    stats.k = numel(Ms);                                                                    % number of groups 
    stats.Ws = 1 ./ Vs;                                                                     % weights = 1 / variances
    
    stats.FE_SumOfWeights = sum(stats.Ws);                                                  % 1) FE sum of weights
    stats.FE_SumOfWeightedMeans = sum(stats.Ws.*Ms);                                        % 2) FE sum of weight x mean
    stats.RE_SumOfEffectxVariance = sum(Ms.^2 ./ Vs);                                       % 9) sum of effect^2 / variance
    stats.FE_SumOfSquaredWeights = sum(stats.Ws.^2);                                        % 10) sum of FE weight^2
    stats.Q = stats.RE_SumOfEffectxVariance - (stats.FE_SumOfWeightedMeans.^2 ./stats.FE_SumOfWeights);% 11) Q = heterogeneity; Cochran's Q-test, Hedges and Olkin 1985
    stats.Qdf = stats.k - 1;                                                                % 12) df = groups - 1
    stats.Qp = 1 - chi2cdf(stats.Q, stats.Qdf);                                             % 13) Chi^2 p-value, of Q with df
    
   
    %% METHODS TO ESTIMATE THE BETWEEN-STUDY VARIANCE, TAU-SQUARED_______
    % CI tau^2 obtained with method of Viechtbauer (2007a)
    % https://doi.org/10.1002/jrsm.1164
    % Tau^2 is unreliable with small N, so CIs are useful
    
    switch method
        case {'FE'}                                                                         % no between-study variance
            %% FIXED EFFECTS META-ANALYSIS_______________________________
            stats.FE_Mean = stats.FE_SumOfWeightedMeans ./ stats.FE_SumOfWeights;           % 3) FE meta (weighted) mean = sum of (weight x M) / sum of (weights)
            stats.FE_SE = sqrt(1./stats.FE_SumOfWeights);                                   % 4) FE meta SE = sqrt(1/(sum of weights))
            stats.FE_CI(1) =  stats.FE_Mean + (stats.FE_SE .* norminv(alpha./2));           % 5) FE meta CI = M + (SE * norminv(alpha))
            stats.FE_CI(2) =  stats.FE_Mean - (stats.FE_SE .* norminv(alpha./2));           % 6) FE meta CI = M - (SE * norminv(alpha))
            stats.FE_Z = (stats.FE_Mean - comparison) ./ stats.FE_SE;                       % 7) FE meta Z-score = M ./ SE
            stats.FE_p = (1-normcdf(abs(stats.FE_Z))).*2;                                   % 8) FE meta p-value from abs Z-score
	    
	    
        %% Categorisation according to Veroniki et al (2015)_____________
	    
	    %% Method of moments estimators__________________________________
	    
	    % Non-iterative____________________________________________________
    
	    case {'DL'}
            % DerSimonian and Laird 1986; Raudenbush 2009__________________
	        % DerSimonian and Laird 1986 (https://doi.org/10.1016/0197-2456(86)90046-2); Raudenbush 2009
	        % criticised for under-estimating true between-study variance, with too-narrow CIs, especially when Tau-squared is large and/or when using 2x2 designs / odds ratios/risk ratios
	        % because negative Tau-squared has to be truncated, the estimator is biased upwards
	        % inefficient (variable) when study sizes vary a lot and Tau-squared is high
	        % can be efficient when k studies is large
            stats.RE_TauSquared = max(0,(stats.Q - stats.Qdf) ./ (stats.FE_SumOfWeights - (stats.FE_SumOfSquaredWeights ./ stats.FE_SumOfWeights))); % 14) Tau^2 = between-study variance
	        stats.RE_ISquared = 100 .* max(0, (stats.Q - stats.Qdf) ./ stats.Q);            % 23) excess variance as percentage of total variance observed in the effect size estimates across studies - akin to intraclass correlation
	        stats.RE_HSquared = stats.Q ./ stats.Qdf;                                       % 24) total variability / within-study variance
    
        case {'DLp'}
	        % Tau-squared can not be zero -> add an arbitrary non-zero positive number (e.g. 0.01) when Tau-squared is zero or negative
	        % how to choose that number?
	    
	    case {'DL2'}
	        % Calculate Q as in DL, then adjust...
	    
	    case {'HO'}
	        % Hedges and Olkin 1985 (https://doi.org/10.1016/B978-0-08-057065-5.50014-2); Cochran; Raudenbush 2009
	        % also variance component type estimator
	        % set sample variance to its expected value, with unweighted mean of y
	        % simple, but not widely used; gives larger estimates than DL and ML and REML
	        % unbiased before truncation, as long as the sampling variances are unbiased
	        % good with large between-study variance
	        % Sy^2 = (1/k-1).*Sum of (yi-ybar)^2
	        % Tau^2 = max(0, Sy^2 - (1/k).*sum(Vi))
	        stats.RE_TauSquared = max(0,((1./stats.Qdf) .* sum((Ms - mean(Ms)).^2)) - ((1./stats.k).*sum(Vs)));
	        stats.RE_ISquared = [];%100 .* max(0, (stats.Q - stats.Qdf) ./ stats.Q);        % 23) excess variance as percentage of total variance observed in the effect size estimates across studies - akin to intraclass correlation
	        stats.RE_HSquared = [];%stats.Q ./ stats.Qdf;                                   % 24) total variability / within-study variance
	        
	    case {'HO2'}
	        % DerSimonian and Kacker (2007; https://doi.org/10.1016/j.cct.2006.04.004)
	        % same as two-step DL, replacing the DL with the HO
	
	    case {'HM'}
	        % Tau-squared can not be zero -> uses quadratic form
	        % for small differences, produces large Tau^2
	        % Tau-squared = Q^2 / ((2*(k-1) + Q) * (sum of FE weights * (sum of squared FE weights / sum of FE weights))
	        stats.RE_TauSquared = stats.Q.^2 ./ ((2.*(stats.Qdf) + stats.Q) .* (stats.FE_SumOfWeights * (stats.FE_SumOfSquaredWeights ./ stats.FE_SumOfWeights))); % THIS HAS NOT BEEN CHECKED..!
    
	    case {'HS'}
	        % Hunter and Schmidt (2004; "Methods of Meta-Analysis")
	        % Tau-squared = max(0, (Q - k) ./ Sum Of FE Weights
	        stats.RE_TauSquared = max(0, (stats.Q - stats.k) ./ stats.FE_SumOfWeights);
	        stats.RE_ISquared = [];%100 .* max(0, (stats.Q - stats.Qdf) ./ stats.Q);        % 23) excess variance as percentage of total variance observed in the effect size estimates across studies - akin to intraclass correlation
	        stats.RE_HSquared = [];%stats.Q ./ stats.Qdf;                                   % 24) total variability / within-study variance
	        
	    % Iterative________________________________________________________
	    case {'PM'}                                                                         % see https://methods.cochrane.org/sites/methods.cochrane.org.statistics/files/uploads/SMG_training_course_2016/cochrane_smg_training_2016_viechtbauer.pdf
	        % equivalent to the empirical Bayes estimator (Morris 1983)
    
    
	    %% Maximum likelihood estimators_________________________________
            
	    % Iterative, non-negative__________________________________________
	    case {'ML'}
	        % Viechtbauer 2005; Raudenbush 2009
	        % https://www.datacamp.com/tutorial/maximum-likelihood-estimation-mle
	        % https://doi.org/10.23943/princeton/9780691137285.003.0010
	        % https://journals.sagepub.com/doi/abs/10.1177/1368430214558311
	        % When the estimate of Tau-squared is also of interest, ML has been recommended
	        % ML depends on the maximisation method, which may fail, especially when k is small
	        % likelihood methods are asymtotically unbiased
	        % ML has a small MSE, but large negative bias for large Tau^2 when k is small and small studies are included
	        % the known properties are under assumption of normal distribution -> non-normal may behave differently
	        % Viechbauer (2005) recommends avoiding the ML
	        % 
	        % log-likelihood function, lnL(mu, Tau^2) = -(k/2)ln(2*pi) - (0.5)*sum(ln(vi + Tau^2) - (0.5)*(sum(yi - mu)^2) / (vi + Tau^2)
	        % 
	        % iteration over mu and Tau^2 to give:
	        % mu = Sum of RE-weighted means / sum of RE weights
	        % Tau^2 = max(0, [sum of squared RE weights * ((yi - mu)^2 - vi) / sum of squared weights])
	        % 
	        % iteration via Newton-Raphson, scoring, simplex, or expectation-maximisation methods
	    
	    case {'REML'}
	        % Viechtbauer 2005; Raudenbush 2009
	        % JASP: "When the focus of the analysis is directed to the fixed effects part of the model, Restricted ML has been recommended."
	        % REML estimator is approximately unbiased and quite efficient
	        % REML can correct for the negative bias of the ML method
	        % REML is preferable when large studies are included
	        % mu is removed from the iteration - one-dimensional search
	        
	        % log-likelihood function, lnL(Tau^2) = -(k/2)ln(2*pi) - (0.5)*sum(ln(vi + Tau^2) - (0.5)*(sum(yi - mu)^2) / (vi + Tau^2) - (0.5)*ln(sum of 1/(vi + Tau^2))
	        % 
	        % iteration over mu and Tau^2 to give:
	        % mu = Sum of RE-weighted means / sum of RE weights
	        % Tau^2 = max(0, [sum of squared RE weights * ((yi - mu)^2 - vi) / sum of squared RE weights + (1/sum of RE weights)])
	        
	    case {'AREML'}
	        % Tau^2 = max(0, [sum of squared RE weights * ( (k/(k-1))*(yi - mu)^2 - vi) / sum of squared weights])
    
	    %%  Model error variance estimators (non-iterative)______________
	    case {'SJ'}                                                         % Sidik and Jonkman 2005a,b
	        % Tau0^2 = sum(yi - ybar)^2 / k
	        % r-hat = vi/Tau0^2 (Tau0^2>0)	for all means
	        % q-hat = r-hat + 1			for all means
	        % mu-hat = sum((1/q-hat).*yi / sum(1./q-hat)
	        % Tau^2 = (1/(k - 1)) * sum((1/q-hat)*(yi - mu-hat)^2
	        Tau0 = sum((Ms - mean(Ms)).^2) ./ stats.k;
	        rhat = Vs./Tau0;
	        qhat = rhat + 1;
	        muhat = sum((1./qhat).*Ms) ./ sum(1./qhat);
	        stats.RE_TauSquared = (1./stats.Qdf) .* sum((1./qhat).*(Ms - muhat).^2);
	        stats.RE_ISquared = [];%100 .* max(0, (stats.Q - stats.Qdf) ./ stats.Q);        % 23) excess variance as percentage of total variance observed in the effect size estimates across studies - akin to intraclass correlation
	        stats.RE_HSquared = [];%stats.Q ./ stats.Qdf;                                   % 24) total variability / within-study variance
	        
	    %% Bayes estimators______________________________________________
	    case {'RB'}
	    
	    case {'RBp'}                                                                        % Tau-squared can not be zero
	    
	    case {'FB'}
	    
	    case {'BM'}
	    
	    case {'EB'}                                                                         % Morris 1983; Berkey et al. 1995; see https://methods.cochrane.org/sites/methods.cochrane.org.statistics/files/uploads/SMG_training_course_2016/cochrane_smg_training_2016_viechtbauer.pdf

	    %% Bootstrap estimator___________________________________________
        case {'Non-parametric bootstrap DerSimonian-Laird','DLb'}
             % Kontopantelis et al. (2013; https://doi.org/10.1371/journal.pone.0069930)
	     % randomly sample B sets of studies with replacement, calculate DL stat, then take mean of DL stats
	     % may only work well with large k, and may be biased with small k
	     % Bootstrap could be applied for all estimators...


        %% OTHERS?_______________________________________________________
	    case {'Q-profile'}
	    
	    case {'Knapp and Hartung'}                                                          % see https://methods.cochrane.org/sites/methods.cochrane.org.statistics/files/uploads/SMG_training_course_2016/cochrane_smg_training_2016_viechtbauer.pdf
	    
	    case {'Mantel-Haenszel'}	                                                        % when inverse-variance weighting is not appropriate, eg in dichotomous outcomes with rare events
	    
	    case {'Peto'}                                                                       % for balanced group sizes
	    
    end
    
    switch method
        case {'DL', 'HO', 'HS', 'SJ'}
            stats.RE_SumOfWeights = sum(1 ./ (Vs + stats.RE_TauSquared));                   % 15) sum of RE weights = 1 / (within variance + between variance)
            stats.RE_SumOfEffectsOverVariances = sum(Ms ./ (Vs + stats.RE_TauSquared));     % 16) sum of Estimate / RE variance (within variance + between variance)
            stats.RE_Mean = stats.RE_SumOfEffectsOverVariances ./ stats.RE_SumOfWeights;    % 17) RE meta-analytic mean
            stats.RE_SE = sqrt(1 ./ stats.RE_SumOfWeights);                                 % 18) RE meta-analytic SE
            stats.RE_CI(1) = stats.RE_Mean + (stats.RE_SE .* norminv(alpha./2));            % 19) RE upper CI
            stats.RE_CI(2) = stats.RE_Mean - (stats.RE_SE .* norminv(alpha./2));            % 20) RE lower CI
            stats.RE_Z = (stats.RE_Mean - comparison) ./ stats.RE_SE;                       % 21) RE Z score
            stats.RE_p = (1-normcdf(abs(stats.RE_Z))).*2;                                   % 22) RE meta p-value from abs Z-score
        
	    otherwise
    end
    
    %% FUNNEL PLOT ASYMMETRY & DIAGNOSTICS_______________________________
    
    % trim and fill________________________________________________________
    
    % rank test (Kendall's Tau)____________________________________________
    
    % regression test (Egger's sei)________________________________________
    
    % Fail-safe N (Rosenthal's file drawer analysis)_______________________
    
end