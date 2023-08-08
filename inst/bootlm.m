% -- Function File: bootlm (Y, GROUP)
% -- Function File: bootlm (Y, GROUP, ..., NAME, VALUE)
% -- Function File: bootlm (Y, GROUP, ..., 'dim', DIM)
% -- Function File: bootlm (Y, GROUP, ..., 'continuous', CONTINUOUS)
% -- Function File: bootlm (Y, GROUP, ..., 'model', MODELTYPE)
% -- Function File: bootlm (Y, GROUP, ..., 'varnames', VARNAMES)
% -- Function File: bootlm (Y, GROUP, ..., 'method', METHOD)
% -- Function File: bootlm (Y, GROUP, ..., 'method', 'bayesian', 'prior', PRIOR)
% -- Function File: bootlm (Y, GROUP, ..., 'alpha', ALPHA)
% -- Function File: bootlm (Y, GROUP, ..., 'display', DISPOPT)
% -- Function File: bootlm (Y, GROUP, ..., 'contrasts', CONTRASTS)
% -- Function File: bootlm (Y, GROUP, ..., 'nboot', NBOOT)
% -- Function File: bootlm (Y, GROUP, ..., 'clustid', CLUSTID)
% -- Function File: bootlm (Y, GROUP, ..., 'blocksz', BLOCKSZ)
% -- Function File: bootlm (Y, GROUP, ..., 'posthoc', POSTHOC)
% -- Function File: bootlm (Y, GROUP, ..., 'seed', SEED)
% -- Function File: STATS = bootlm (...)
% -- Function File: [STATS, BOOTSTAT] = bootlm (...)
% -- Function File: [STATS, BOOTSTAT, AOVSTAT] = bootlm (...)
% -- Function File: [STATS, BOOTSTAT, AOVSTAT, X] = bootlm (...)
% -- Function File: [STATS, BOOTSTAT, AOVSTAT, X, L] = bootlm (...)
%
%        Fits a linear model with categorical and/or continuous predictors (i.e.
%     independent variables) on a continuous outcome (i.e. dependent variable)
%     and computes the following statistics for each regression coefficient:
%          - name: the name(s) of the regression coefficient(s)
%          - coeff: the value of the regression coefficient(s)
%          - CI_lower: lower bound(s) of the 95% confidence interval (CI)
%          - CI_upper: upper bound(s) of the 95% confidence interval (CI)
%          - p-val: two-tailed p-value(s) for the parameter(s) being equal to 0
%        By default, confidence intervals and Null Hypothesis Significance Tests
%     (NHSTs) for the regression coefficients (H0 = 0) are calculated by wild
%     bootstrap-t and so are robust when normality and homoscedasticity cannot
%     be assumed.
%
%        Usage of this function is very similar to that of 'anovan'. Data (Y)
%     is a single vector y with groups specified by a corresponding matrix or
%     cell array of group labels GROUP, where each column of GROUP has the same
%     number of rows as Y. For example, if 'Y = [1.1;1.2]; GROUP = [1,2,1; 
%     1,5,2];' then observation 1.1 was measured under conditions 1,2,1 and
%     observation 1.2 was measured under conditions 1,5,2. If the GROUP provided
%     is empty, then the linear model is fit with just the intercept (i.e. no
%     predictors).
%
%     'bootlm' can take a number of optional parameters as name-value
%     pairs.
%
%     '[...] = bootlm (Y, GROUP, ..., 'dim', DIM)'
%
%       <> DIM is a scalar or vector specifying the dimension(s) over which
%          'bootlm' calculates and returns estimated marginal means instead of
%          regression coefficients. For example, the value [1 3] computes the
%          estimated marginal mean for each combination of the levels of the
%          first and third predictors. The default value is empty, which makes
%          'bootlm' return the statistics for for the model coefficients. If DIM
%          is, or includes, a continuous predictor then 'bootlm' will return an
%          error. The following statistics are printed when specifying 'dim':
%             - name: the name(s) of the estimated marginal mean(s)
%             - mean: the estimated marginal mean(s)
%             - CI_lower: lower bound(s) of the 95% confidence interval (CI)
%             - CI_upper: upper bound(s) of the 95% confidence interval (CI)
%             - n: the sample size used to estimate the mean
%
%     '[...] = bootlm (Y, GROUP, ..., 'continuous', CONTINUOUS)'
%
%       <> CONTINUOUS is a vector of indices indicating which of the
%          columns (i.e. predictors) in GROUP should be treated as
%          continuous predictors rather than as categorical predictors.
%          The relationship between continuous predictors and the outcome
%          should be linear.
%
%     '[...] = bootlm (Y, GROUP, ..., 'model', MODELTYPE)'
%
%       <> MODELTYPE can specified as one of the following:
%
%             o 'linear' (default): compute N main effects with no
%               interactions.
%
%             o 'interaction': compute N effects and N*(N-1) interactions
%
%             o 'full': compute the N main effects and interactions at
%               all levels
%
%             o a scalar integer: representing the maximum interaction
%               order
%
%             o a matrix of term definitions: each row is a term and
%               each column is a predictor
%
%               -- Example:
%               A two-way design with interaction would be: [1 0; 0 1; 1 1]
%
%     '[...] = bootlm (Y, GROUP, ..., 'varnames', VARNAMES)'
%
%       <> VARNAMES must be a cell array of strings with each element
%          containing a predictor name for each column of GROUP. By default
%          (if not parsed as optional argument), VARNAMES are
%          'X1','X2','X3', etc.
%
%     '[...] = bootlm (Y, GROUP, ..., 'method', METHOD)'
%
%       <> METHOD can be specified as one of the following:
%
%             o 'wild' (default): Wild bootstrap-t, using the 'bootwild'
%               function. Please see the help documentation for the function
%               'bootwild' for more information about this method.
%
%             o 'bayesian': Bayesian bootstrap, using the 'bootbayes' function.
%                Please see the help documentation below in the function
%               'bootbayes' for more information about this method.
%
%             Note that p-values are a frequentist concept and are only computed
%             and returned from bootlm when the METHOD is 'wild'.
%
%     '[...] = bootlm (Y, GROUP, ..., 'method', 'bayesian', 'prior', PRIOR)'
%
%       <> Sets the prior for Bayesian bootstrap. Possible values are:
%
%             o scalar: A positive real numeric scalar to parametrize
%                  the form of the symmetric Dirichlet distribution. The
%                  Dirichlet distribution is the conjugate PRIOR used to
%                  randomly generate weights for linear least squares fitting
%                  of the observed data and subsequently to estimate the
%                  posterior for the regression coefficients by nonparametric
%                  Bayesian bootstrap.
%
%             o 'auto': Sets a value for PRIOR that effectively incorporates
%                  Bessel's correction a priori such that the variance of the
%                  posterior (i.e. the rows of BOOTSTAT) becomes an unbiased
%                  estimator of the sampling variance. The calculation used for
%                  'auto' is as follows:
% 
%                     PRIOR = 1 - 2 / N
% 
%                  For block or cluster bootstrap, N corresponds to the number
%                  of blocks or clusters (i.e. the number of independent
%                  sampling units).
%                       The 'auto' setting is recommended but is only available
%                  for Bayesian bootstrap of the estimated marginal means and
%                  for the posthoc tests (not the regression coefficients).
%
%               The default value of PRIOR is the scalar: 1, which corresponds
%               to Bayes rule: a uniform (or flat) Dirichlet distribution
%               (over all points in its support). Please see the help
%               documentation for the function 'bootbayes' for more information
%               about the prior.
%
%     '[...] = bootlm (Y, GROUP, ..., 'alpha', ALPHA)'
%
%       <> ALPHA is numeric and sets the lower and upper bounds of the
%          confidence or credible interval(s). The value(s) of ALPHA must be
%          between 0 and 1. ALPHA can either be:
%
%             o scalar: Set the central mass of the intervals to 100*(1-ALPHA)%.
%                  For example, 0.05 for a 95% interval. If METHOD is 'wild',
%                  then the intervals are symmetric bootstrap-t confidence
%                  intervals. If METHOD is 'bayesian', then the intervals are
%                  shortest probability credible intervals.
%
%             o vector: A pair of probabilities defining the lower and upper
%                  and upper bounds of the interval(s) as 100*(ALPHA(1))% and 
%                  100*(ALPHA(2))% respectively. For example, [.025, .975] for
%                  a 95% interval. If METHOD is 'wild', then the intervals are
%                  assymmetric bootstrap-t confidence intervals. If METHOD is
%                  'bayesian', then the intervals are simple percentile credible
%                  intervals.
%
%               The default value of ALPHA is the scalar: 0.05.
%
%     '[...] = bootlm (Y, GROUP, ..., 'display', DISPOPT)'
%
%       <> DISPOPT can be either 'on' (or true, default) or 'off' (or false)
%          and controls the display of the model formula, a table of model
%          parameter estimates and a figure of diagnostic plots. The p-values
%          are formatted in APA-style.
%
%     '[...] = bootlm (Y, GROUP, ..., 'contrasts', CONTRASTS)'
%
%       <> CONTRASTS can be specified as one of the following:
%
%             o A string corresponding to one of the built-in contrasts
%               listed below:
%
%                  o 'anova' or 'simple' (default): Simple (ANOVA) contrast
%                    coding. The intercept represents the grand mean. Each 
%                    slope coefficient represents the difference between one
%                    level of a predictor (or interaction between predictors) to
%                    the first level for that/those predictor(s), averaged over
%                    all levels of the other predictor(s). The first (or
%                    reference level) of the predictor(s) is defined as the
%                    first level of the predictor (or combination of the
%                    predictors) listed in the GROUP argument. The columns of
%                    this contrast coding scheme sum to zero. This type of
%                    contrast is ideal for nominal predictor variables that
%                    have an obvious reference or control group and that are
%                    modelled together with a covariate or blocking factor.
%
%                  o 'poly': Polynomial contrast coding for trend analysis.
%                    The intercept represents the grand mean. The remaining
%                    slope coefficients returned are for linear, quadratic,
%                    cubic etc. trends across the levels. In addition to the
%                    columns of this contrast coding scheme summing to zero,
%                    this contrast coding is orthogonal (i.e. the off-diagonal
%                    elements of its autocovariance matrix are zero) and so
%                    the slope coefficients are independent. This type of
%                    contrast is ideal for ordinal predictor variables, in
%                    particular, predictors with ordered levels that are evenly
%                    spaced.
%
%                  o 'helmert': Helmert contrast coding. The intercept
%                    represents the grand mean. Each slope coefficient
%                    represents the difference between one level of a predictor
%                    (or interaction between predictors) with the mean of the
%                    subsequent levels, where the order of the predictor levels
%                    is as they appear in the GROUP argument. In addition to the
%                    columns of this contrast coding scheme summing to zero,
%                    this contrast coding is orthogonal (i.e. the off-diagonal
%                    elements of its autocovariance matrix are zero) and so the
%                    slope coefficients are independent. This type of contrast
%                    is ideal for predictor variables that are either ordinal,
%                    or nominal with their levels ordered such that the contrast
%                    coding reflects tests of some hypotheses of interest about
%                    the nested grouping of the predictor levels.
%
%                  o 'effect': Deviation effect coding. The intercept represents
%                    the grand mean. Each slope coefficient compares one level
%                    of a predictor (or interaction between predictors) with the
%                    grand mean. Note that a slope coefficient is omitted for
%                    the first level of the predictor(s) listed in the GROUP
%                    argument. The columns of this contrast coding scheme sum to
%                    zero. This type of contrast is ideal for nominal predictor
%                    variables when there is no obvious reference group.
%
%                  o 'sdif' or 'sdiff': Successive differences contrast coding.
%                    The intercept represents the grand mean. Each slope
%                    coefficient represents the difference between one level of
%                    a predictor (or interaction between predictors) to the
%                    previous one, where the order of the predictor levels is
%                    as they appear in the GROUP argument. The columns of this
%                    contrast coding coding scheme sum to zero. This type of
%                    contrast is ideal for ordinal predictor variables.
%
%                  o 'treatment': Treatment contrast (or dummy) coding. The
%                    intercept represents the mean of the first level of all
%                    the predictors. Each slope coefficient compares one
%                    level of a predictor (or interaction between predictors)
%                    with the first level for that/those predictor(s), at the
%                    first level of all the other predictors. The first (or
%                    reference level) of the predictor(s) is defined as the
%                    first level of the predictor (or combination of the
%                    predictors) listed in the GROUP argument. This type of
%                    contrast is ideal for one-way designs or factorial designs
%                    of nominal predictor variables that have an obvious
%                    reference or control group.
%
%            <> A matrix containing a custom contrast coding scheme (i.e.
%               the generalized inverse of contrast weights). Rows in
%               the contrast matrices correspond to predictor levels in the
%               order that they first appear in the GROUP column. The
%               matrix must contain the same number of columns as there
%               are the number of predictor levels minus one.
%
%          If the linear model contains more than one predictor and a
%          built-in contrast coding scheme was specified, then those
%          contrasts are applied to all predictors. To specify different
%          contrasts for different predictors in the model, CONTRASTS should
%          be a cell array with the same number of cells as there are
%          columns in GROUP. Each cell should define contrasts for the
%          respective column in GROUP by one of the methods described
%          above. If cells are left empty, then the default contrasts
%          are applied. Contrasts for cells corresponding to continuous
%          predictors are ignored.
%
%     '[...] = bootlm (Y, GROUP, ..., 'nboot', NBOOT)'
%
%       <> Specifies the number of bootstrap resamples, where NBOOT must be a
%          positive integer. If empty, the default value of NBOOT is 9999.
%
%     '[...] = bootlm (Y, GROUP, ..., 'clustid', CLUSTID)'
%
%       <> Specifies a vector or cell array of numbers or strings respectively
%          to be used as cluster labels or identifiers. Rows of the data with
%          the same CLUSTID value are treated as clusters with dependent errors.
%          If empty (default), no clustered resampling is performed and all
%          errors are treated as independent. The standard errors computed are
%          cluster robust.
%
%     '[...] = bootlm (Y, GROUP, ..., 'blocksz', BLOCKSZ)'
%
%       <> Specifies a scalar, which sets the block size for bootstrapping when
%          the errors have serial dependence. Rows of the data within the same
%          block are treated as having dependent errors. If empty (default),
%          no clustered resampling is performed and all errors are treated
%          as independent. The standard errors computed are cluster robust.
%
%     '[...] = bootlm (Y, GROUP, ..., 'posthoc', POSTHOC)'
%
%       <> When DIM is specified, POSTHOC comparisons along DIM can be one of
%          the following:
%
%             o 'none' (default): No posthoc comparisons are performed. The
%               statistics returned are for the estimated marginal means.
%
%             o 'pairwise' : Pairwise comparisons are performed.
%
%             o 'trt_vs_ctrl' : Treatment vs. Control comparisons are performed.
%                The control is the first group number 1 as returned when 
%                POSTHOC is set to 'none'.
%
%             o {'trt_vs_ctrl', k} : Treatment vs. Control comparisons are
%                performed. The control is group number k as returned when
%                POSTHOC is set to 'none'.
%
%          All of the posthoc comparisons use the Holm-Bonferroni procedure to
%          control the type I error rate. The confidence intervals are not
%          adjusted for multiple comparisons.
%
%     '[...] = bootlm (Y, GROUP, ..., 'seed', SEED)' initialises the Mersenne
%     Twister random number generator using an integer SEED value so that
%     'bootlm' results are reproducible.
%
%     'bootlm' can return up to four output arguments:
%
%     'STATS = bootlm (...)' returns a structure with the following fields:
%        - 'method': The bootstrap method
%        - 'name': The names of each of the estimates
%        - 'estimate': The value of the estimates
%        - 'CI_lower': The lower bound(s) of the confidence/credible interval(s)
%        - 'CI_upper': The upper bound(s) of the confidence/credible interval(s)
%        - 'pval': The p-value(s) for the hypothesis that the estimate(s) = 0
%        - 'fpr': The false positive risk 
%        - 'n': The sample size(s)
%        - 'prior': The prior used for Bayesian bootstrap
%
%        Note that the p-values returned are truncated at the resolution
%        limit determined by the number of bootstrap replicates, specifically 
%        1 / (NBOOT + 1).
%
%     '[STATS, BOOTSTAT] = bootlm (...)' also returns a P x NBOOT matrix of
%     bootstrap statistics for the estimated parameters, where P is the number
%     of parameters estimated in the model. Depending on the DIM and POSTHOC
%     input arguments set by the user, the estimated parameters whose bootstrap
%     statistics are returned will be either regression coefficients, the
%     estimated marginal means, or the mean differences between groups of a
%     categorical predictor for posthoc testing.
%
%     '[STATS, BOOTSTAT, AOVSTAT] = bootlm (...)' also computes and returns
%     bootstrapped ANOVA statistics in a structure with the following fields: 
%        - 'MODEL': The formula of the linear model(s) in Wilkinson's notation
%        - 'SS': Sum-of-squares
%        - 'DF': Degrees of freedom
%        - 'MS': Mean-squares
%        - 'F': F-Statistic
%        - 'PVAL': p-values
%        - 'SSE': Sum-of-Squared Error
%        - 'DFE': Degrees of Freedom for Error
%        - 'MSE': Mean Squared Error
%     The ANOVA implemented uses sequential (type I) sums-of-squares and so the
%     results and their interpretation depend on the order of predictors in the
%     GROUP variable (when the design is not balanced). Thus, the null model
%     used for comparison for each model is the model listed directly above it
%     in AOVSTAT; for the first model, the null model is the intercept-only
%     model. Note that ANOVA statistics are only returned for wild bootstrap
%     AND when no other statistics are requested (i.e. estimated marginal means
%     or posthoc tests). The bootstrap is achieved by Wild bootstrap of the
%     residuals from the full model.
%
%     '[STATS, BOOTSTAT, AOVSTAT, X] = bootlm (...)' also returns the design
%     matrix for the linear  model.
%
%     '[STATS, BOOTSTAT, AOVSTAT, X, L] = bootlm (...)' also returns the
%     hypothesis matrix used to compute the estimated marginal means or posthoc
%     tests from the regression coefficients.
%
%  bootlm (version 2023.08.02)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/
%
%  Copyright 2019 Andrew Charles Penn
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%  GNU General Public License for more details.

function [STATS, BOOTSTAT, AOVSTAT, X, L] = bootlm (Y, GROUP, varargin)

    if (nargin < 2)
      error (cat (2, 'bootlm usage: ''bootlm (Y, GROUP)''; ', ...
                     ' atleast 2 input arguments required'))
    end
    if (nargout > 5)
      error ('bootlm: Too many output arguments')
    end

    % Check if running in Octave (else assume Matlab)
    info = ver; 
    ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

    % Check supplied parameters
    if ((numel (varargin) / 2) ~= fix (numel (varargin) / 2))
      error ('bootlm: wrong number of arguments.')
    end
    MODELTYPE = 'linear';
    DISPLAY = 'on';
    VARNAMES = [];
    CONTINUOUS = [];
    CONTRASTS = {};
    ALPHA = 0.05;
    DIM = [];
    NBOOT = 9999;
    SEED = rand ('seed');
    DEP = [];
    POSTHOC = 'none';
    METHOD = 'wild';
    PRIOR = 1;
    L = [];
    AOVSTAT = [];
    for idx = 3:2:nargin
      name = varargin{idx-2};
      value = varargin{idx-1};
      switch (lower (name))
        case 'model'
          MODELTYPE = value;
        case 'continuous'
          CONTINUOUS = value;
        case 'random'
          warning (cat (2, 'bootlm: ''random'' name-value pair is not', ...
                           ' supported and will be ignored.'))
        case 'nested'
          error (cat (2, 'bootlm: ''nested'' name-value pair is not', ...
                         ' supported. Please use ''CLUSTID'' or', ...
                         ' ''BLOCKSZ'' input argument.'))
        case 'sstype'
          warning (cat (2, 'bootlm: ''sstype'' name-value pair is not', ...
                           ' supported and will be ignored.'))
        case 'varnames'
          VARNAMES = value;
        case 'display'
          DISPLAY = value;
        case 'contrasts'
          CONTRASTS = value;
        case 'alpha'
          ALPHA = value;
        case {'clustid', 'blocksz'}
          DEP = value;
        case {'dim', 'dimension'}
          DIM = value;
        case {'posthoc', 'posttest'}
          POSTHOC = value;
        case 'nboot'
          NBOOT = value;
        case 'method'
          METHOD = value;
        case 'prior'
          PRIOR = value;
        case 'seed'
          SEED = value;
        otherwise
          error (sprintf ('bootlm: parameter %s is not supported', name))
      end
    end

    % Most error checking for NBOOT, ALPHA and SEED is handled by the functions
    % bootwild and bootbayes
    if (size (ALPHA,1) > 1)
      ALPHA = ALPHA.';
    end

    % Evaluate continuous input argument
    if (isnumeric (CONTINUOUS))
      if (any (CONTINUOUS ~= abs (fix (CONTINUOUS))))
        error (cat (2, 'bootlm: the value provided for the CONTINUOUS', ...
                       ' parameter must be a positive integer'))
      end
    else
      error (cat (2, 'bootlm: the value provided for the CONTINUOUS', ...
                     ' parameter must be numeric'))
    end

    % Accomodate for different formats for GROUP
    % GROUP can be a matrix of numeric identifiers of a cell arrays
    % of strings or numeric identifiers
    N = size (GROUP, 2); % number of predictors
    n = numel (Y);       % total number of observations
    if (prod (size (Y)) ~= n)
      error ('bootlm: for ''bootlm (Y, GROUP)'', Y must be a vector')
    end
    if (numel (unique (CONTINUOUS)) > N)
      error (cat (2, 'bootlm: the number of predictors assigned as', ...
                     ' continuous cannot exceed the number of', ...
                     ' predictors in GROUP'))
    end
    if (any ((CONTINUOUS > N) | any (CONTINUOUS <= 0)))
      error (cat (2, 'bootlm: one or more indices provided in the value', ...
                     ' for the continuous parameter are out of range'))
    end
    cont_vec = false (1, N);
    cont_vec(CONTINUOUS) = true;
    if (iscell (GROUP))
      if (size (GROUP, 1) == 1)
        tmp = cell (n, N);
        for j = 1:N
          if (isnumeric (GROUP{j}))
            if (ismember (j, CONTINUOUS))
              tmp(:,j) = num2cell (GROUP{j});
            else
              tmp(:,j) = cellstr (num2str (GROUP{j}));
            end
          else
            if (ismember (j, CONTINUOUS))
              error ('bootlm: continuous predictors must be a numeric datatype')
            end
            tmp(:,j) = GROUP{j};
          end
        end
        GROUP = tmp;
      end
    end
    if (~ isempty (GROUP))
      if (size (GROUP,1) ~= n)
        error (cat (2, 'bootlm: GROUP must be a matrix with the same', ...
                       ' number of rows as Y'))
      end
    end
    if (~ isempty (VARNAMES))
      if (iscell (VARNAMES))
        if (all (cellfun (@ischar, VARNAMES)))
          nvarnames = numel(VARNAMES);
        else
          error (cat (2, 'bootlm: all variable names must be character', ...
                         ' or character arrays'))
        end
      elseif (ischar (VARNAMES))
        nvarnames = 1;
        VARNAMES = {VARNAMES};
      elseif (isstring (VARNAMES))
        nvarnames = 1;
        VARNAMES = {char(VARNAMES)};
      else
        error (cat (2, 'bootlm: varnames is not of a valid type. Must be', ...
               ' a cell array of character arrays, character array or string'))
      end
    else
      nvarnames = N;
      VARNAMES = arrayfun(@(x) ['X',num2str(x)], 1:N, 'UniformOutput', 0);
    end
    if (nvarnames ~= N)
      error (cat (2, 'bootlm: number of variable names is not equal', ...
                     ' to the number of grouping variables'))
    end

    % Evaluate contrasts (if applicable)
    if isempty (CONTRASTS)
      CONTRASTS = cell (1, N);
      planned = false;
    else
      if (ischar(CONTRASTS))
        contr_str = CONTRASTS;
        CONTRASTS = cell (1, N);
        CONTRASTS(:) = {contr_str};
      end
      if (~ iscell (CONTRASTS))
        CONTRASTS = {CONTRASTS};
      end
      for i = 1:N
        if (~ isempty (CONTRASTS{i}))
          if (isnumeric(CONTRASTS{i}))
            % Check whether all the columns sum to 0
            if (any (abs (sum (CONTRASTS{i})) > eps ('single')))
              warning (sprintf ( ...
              'Note that the CONTRASTS for predictor %u do not sum to zero', i))
            end
            % Check whether contrasts are orthogonal
            if (any (abs (reshape (corr (CONTRASTS{i}) - ...
                                     eye (size (CONTRASTS{i}, 2)), [], 1)) ...
                                     > eps ('single')))
              warning (sprintf ( ...
              'Note that the CONTRASTS for predictor %u are not orthogonal', i))
            end
          else
            if (~ ismember (lower (CONTRASTS{i}), ...
                            {'simple','anova','poly','helmert','effect', ...
                              'sdif','sdiff','treatment'}))
              error (cat (2, 'bootlm: valid built-in contrasts are:', ...
                            ' ''simple'', ''poly'', ''helmert'',', ...
                            '''effect'', ''sdif'' or ''treatment'''))
            end
          end
        end
      end
      planned = true;
    end

    % Remove NaN or non-finite observations
    if (isempty (GROUP))
      excl = any ([isnan(Y), isinf(Y)], 2);
    else
      XC = GROUP(:,CONTINUOUS);
      if iscell(XC)
        XC = cell2mat (XC);
      end
      excl = any ([isnan(Y), isinf(Y), any(isnan(XC),2), any(isinf(XC),2)], 2);
      GROUP(excl,:) = [];
    end
    Y(excl) = [];
    if (size (Y, 1) == 1)
      Y = Y.';         % If Y is a row vector, make it a column vector
    end
    n = numel (Y);     % Recalculate total number of observations

    % Evaluate model type input argument and create terms matrix if not provided
    msg = cat (2, 'bootlm: the number of columns in the term definitions', ...
                  ' cannot exceed the number of columns of GROUP');
    if (ischar (MODELTYPE))
      switch (lower (MODELTYPE))
        case 'linear'
          MODELTYPE = 1;
        case {'interaction','interactions'}
          MODELTYPE = 2;
        case 'full'
          MODELTYPE = N;
        otherwise
          error ('bootlm: model type not recognised')
      end
    end
    if (isscalar (MODELTYPE))
      TERMS = cell (MODELTYPE,1);
      v = false (1, N);
      switch (lower (MODELTYPE))
        case 1
          % Create term definitions for an additive linear model
          TERMS = eye (N);
        case 2
          % Create term definitions for a model with two predictor interactions
          if (N > 1)
            Nx = nchoosek (N, 2);
          else
            Nx = 0;
          end
          TERMS = zeros (N + Nx, N);
          TERMS(1:N,:) = eye (N);
          for j = 1:N
            for i = j:N-1
              TERMS(N+j+i-1,j) = 1;
              TERMS(N+j+i-1,i+1) = 1;
            end
          end
        otherwise
          if (MODELTYPE > N)
            error (msg);
          end
          % Create term definitions for a full model
          Nx = zeros (1, N-1);
          Nx = 0;
          for k = 1:N
            Nx = Nx + nchoosek(N,k);
          end
          for j = 1:MODELTYPE
            v(1:j) = 1;
            TERMS{j} = flipud (unique (perms (v), 'rows'));
          end
          TERMS = cell2mat (TERMS);
      end
      TERMS = logical (TERMS);
    else
      % Assume that the user provided a suitable matrix of term definitions
      if (size (MODELTYPE, 2) > N)
        error (msg);
      end
      if (~ all (ismember (MODELTYPE(:), [0,1])))
        error (cat (2, 'bootlm: elements of the model terms matrix', ...
                       ' must be either 0 or 1'))
      end
      TERMS = logical (MODELTYPE);
    end
    % Evaluate terms matrix
    Ng = sum (TERMS, 2);
    if (any (diff (Ng) < 0))
      error (cat (2, 'bootlm: the model terms matrix must list main', ...
                     ' effects above/before interactions'))
    end
    % Evaluate terms
    Nm = sum (Ng == 1);
    Nx = sum (Ng > 1);
    Nt = Nm + Nx;
    if (any (any (TERMS(1:Nm,:), 1) ~= any (TERMS, 1)))
      error (cat (2, 'bootlm: all predictors involved in interactions', ...
                     ' must have a main effect'))
    end

    % Create design matrix
    [X, grpnames, nlevels, df, coeffnames, gid, CONTRASTS, ...
     center_continuous] = mDesignMatrix (GROUP, TERMS, ...
     CONTINUOUS, CONTRASTS, VARNAMES, n, Nm, Nx, Ng, cont_vec);
    dft = n - 1;
    dfe = dft - sum (df);
    if (dfe < 1)
      error (cat (2, 'bootlm: there are no error degrees of freedom in', ...
                     ' the specified model'))
    end

    % If applicable, create hypothesis matrix, names and compute sample sizes
    if (isempty (DIM))
      L = 1;
    else
      if (any (DIM < 1))
        error ('bootlm: DIM must contain positive integers')
      end
      if (~ all (ismember (DIM, (1:Nm))))
        error ('bootlm: values in DIM cannot exceed the number of predictors')
      end
      H = X;
      ridx = ~ ismember ((1 : Nm), DIM);
      for i = 1:Nt
        if ( any (and (TERMS(i,:), ridx)) )
          H{i+1}(:,:) = 0;
        end
      end
      H = cell2mat (H);
      L = unique_stable (H, 'rows')';
    end

    % Fit linear model
    X = cell2mat (X);
    [b, sse, resid, ucov, hat] = lmfit (X, Y, ISOCTAVE);

    % Prepare model formula
    TERMNAMES = arrayfun (@(i) sprintf (':%s', VARNAMES{TERMS(i,:)}), ...
                                        (1:Nt), 'UniformOutput', false);
    formula = cell (Nt, 1);
    for i = 1:Nt
      if (i > 1)
        formula{i} = sprintf ('%s + %s', formula{i-1}, TERMNAMES{i}(2:end));
      else 
        formula{1} = sprintf ('Y ~ 1 + %s', TERMNAMES{1}(2:end));
      end
    end

    % Evaluate the dependence structure
    if (isempty (DEP))
      IC = (1:n)';
      IA = IC;
    else
      if (isscalar (DEP))
        % Blocks
        blocksz = DEP;
        G = fix (n / blocksz);
        IC = (G + 1) * ones (n, 1);
        IC(1 : blocksz * G, :) = reshape (ones (blocksz, 1) * (1:G), [], 1);
        [jnk, IA] = unique (IC, 'first');
      else
        % Clusters
        [jnk, IA, IC] = unique_stable (DEP);
        if ( any (size (IC) ~= [n, 1]) )
          error (cat (2, 'bootlm: CLUSTID must be a column vector with the', ... 
                         ' same number of rows as Y'))
        end
      end
    end

    % Use bootstrap methods to calculate statistics
    if isempty (DIM)

      % Error checking
      if ( (~ isempty (POSTHOC)) && (~ strcmpi (lower (POSTHOC), 'none')) )
        error (cat (2, 'bootlm: for posthoc tests you must specify a', ...
                       ' categorical predictor using the DIM input argument'))
      end

      % Model coefficients
      switch (lower (METHOD))
        case 'wild'
          % Perform regression on full model using the specified contrasts
          [STATS, BOOTSTAT] = bootwild (Y, X, DEP, NBOOT, ALPHA, SEED, [], ...
                                        ISOCTAVE);
          % Tidy up
          STATS = rmfield (STATS, {'std_err', 'tstat', 'sse'});
          STATS.n = n;
          STATS.prior = [];
          if (nargout > 2)
            % Perform ANOVA
            AOVSTAT = bootanova (Y, X, cat (1, 1, df), dfe, DEP, NBOOT, ALPHA, ...
                                 SEED, ISOCTAVE);
            AOVSTAT.MODEL = formula;
          end
        case {'bayes', 'bayesian'}
          [STATS, BOOTSTAT] = bootbayes (Y, X, DEP, NBOOT, ...
                                         fliplr (1 - ALPHA), PRIOR, SEED, ...
                                         [], ISOCTAVE);
          % Clean-up
          STATS = rmfield (STATS, {'median', 'bias', 'stdev'});
          STATS.pval = [];
          STATS.fpr = [];
          STATS.n = n;
        otherwise
          error (cat (2, 'bootlm: unrecignised bootstrap method. Use', ...
                         ' ''wild'' or bayesian''.'))
      end

      % Assign names for model coefficients
      NAMES = vertcat (coeffnames{:});
      STATS.name = NAMES;

    else

      % Error checking
      % Check what type of factor is requested in DIM
      if (any (nlevels(DIM) < 2))
          error (cat (2, 'bootlm: DIM must specify only categorical', ...
                         ' factors with 2 or more degrees of freedom.'))
      end
      % Check that all continuous variables were centered
      msg = 'bootlm: model must be refit with a sum-to-zero contrast coding';
      if (any (cont_vec - center_continuous))
        error (msg)
      end
      % Check that the columns of CONTRASTS sum to 0
      for j = 1:N
        if (isnumeric (CONTRASTS{j}))
          if (any (abs (sum (CONTRASTS{j})) > eps('single')))
            error (msg)
          end
        end
      end

      % Create names for estimated marginal means
      idx = cellfun (@(l) find (all (bsxfun (@eq, H, l), 2), 1), ...
                     num2cell (L', 2));
      Np = size (L, 2);
      Nd = numel (DIM);
      NAMES = cell (Np, 1);
      for i = 1 : Np
        str = '';
        for j = 1 : Nd
          str = sprintf('%s%s=%s, ', str, ...
                    num2str (VARNAMES{DIM(j)}), ...
                    num2str (grpnames{DIM(j)}{gid(idx(i),DIM(j))}));
        end
        NAMES{i} = str(1:end-2);
        str = '';
      end

      % Compute sample sizes for each level along dimenion DIM
      U = unique_stable (gid(:,DIM), 'rows');
      n_dim = cellfun (@(u) sum (all (gid(:,DIM) == u, 2)), num2cell (U, 2));

      % Compute number of independent sampling units at each level of DIM
      if (isempty (DEP))
        N_dim = n_dim;
      else
        UC = unique_stable (cat (2, gid(:,DIM), IC), 'rows');
        N_dim = cellfun (@(u) sum (all (UC(:,1:Nd) == u, 2)), num2cell (U, 2));
      end

      switch (lower (POSTHOC))
        case 'none'

          % Model estimated marginal means
          switch (lower (METHOD))
            case 'wild'
              [STATS, BOOTSTAT] = bootwild (Y, X, DEP, NBOOT, ALPHA, SEED, ...
                                            L, ISOCTAVE);
              % Clean-up
              STATS = rmfield (STATS, {'std_err', 'tstat', 'sse'});
              STATS.prior = [];
            case {'bayes', 'bayesian'}
              switch (lower (PRIOR))
                case 'auto'
                  PRIOR = 1 - 2 ./ N_dim;
                  [STATS, BOOTSTAT] = arrayfun (@(i) bootbayes (Y, X, ...
                         DEP, NBOOT, fliplr (1 - ALPHA), PRIOR(i), SEED, ...
                         L(:, i), ISOCTAVE), (1:Np)', 'UniformOutput', false);
                  STATS = flatten_struct (cell2mat (STATS));
                  BOOTSTAT = cell2mat (BOOTSTAT);
                otherwise
                  [STATS, BOOTSTAT] = bootbayes (Y, X, DEP, NBOOT, ...
                                         fliplr (1 - ALPHA), PRIOR, SEED, ...
                                         L, ISOCTAVE);
              end
              % Clean-up
              STATS = rmfield (STATS, {'median', 'bias', 'stdev'});
              STATS.pval = [];
              STATS.fpr = [];
            otherwise
              error (cat (2, 'bootlm: unrecignised bootstrap method. Use', ...
                             ' ''wild'' or bayesian''.'))
          end

          % Add sample sizes to the output structure
          STATS.n = n_dim;

          % Assign NAMES of groups to output structure
          STATS.name = NAMES;

        otherwise

          % Model posthoc comparisons
          if (iscell (POSTHOC))
            if (~ strcmpi (POSTHOC{1}, 'trt_vs_ctrl'))
              error (cat (2, 'bootlm: REF can only be used to specify a', ...\
                             ' control group for ''trt_vs_ctrl'''))
            end
            [L, pairs] = feval (POSTHOC{1}, L, POSTHOC{2:end});
            POSTHOC = POSTHOC{1};
          else
            if (~ ismember (POSTHOC, {'pairwise', 'trt_vs_ctrl'}))
              error (cat (2, 'bootlm: available options for POSTHOC are', ...
                             ' ''pairwise'' and ''trt_vs_ctrl'''))
            end
            [L, pairs] = feval (POSTHOC, L);
          end
          switch (lower (METHOD))
            case 'wild'
              [STATS, BOOTSTAT] = bootwild (Y, X, DEP, NBOOT, ALPHA, SEED, ...
                                            L, ISOCTAVE);
              % Control the type 1 error rate across multiple comparisons
              STATS.pval = holm (STATS.pval);
              % Clean-up
              STATS = rmfield (STATS, {'std_err', 'tstat', 'sse'});
              STATS.prior = [];
            case {'bayes', 'bayesian'}
              switch (lower (PRIOR))
                case 'auto'
                  wgt = bsxfun (@rdivide, N_dim(pairs')', ...
                                sum (N_dim(pairs')', 2));
                  PRIOR = sum ((1 - wgt) .* (1 - 2 ./ N_dim(pairs')'), 2);
                  [STATS, BOOTSTAT] = arrayfun (@(i) bootbayes (Y, X, DEP, ...
                     NBOOT, fliplr (1 - ALPHA), PRIOR(i), SEED, L(:, i), ...
                     ISOCTAVE), (1:size (L, 2))', 'UniformOutput', false);
                  STATS = flatten_struct (cell2mat (STATS));
                  BOOTSTAT = cell2mat (BOOTSTAT);
                otherwise
                  [STATS, BOOTSTAT] = bootbayes (Y, X, DEP, NBOOT, ...
                            fliplr (1 - ALPHA), PRIOR, SEED, L, ISOCTAVE);
              end
              % Clean-up
              STATS = rmfield (STATS, {'median', 'bias', 'stdev'});
              STATS.pval = [];
              STATS.fpr = [];
            otherwise
              error (cat (2, 'bootlm: unrecignised bootstrap method.', ...
                             ' Use ''wild'' or bayesian''.'))
          end

          % Add sample sizes to the output structure
          STATS.n = sum (N_dim(pairs')', 2);

          % Create names of posthoc comparisons and assign to the output
          STATS.name = arrayfun (@(i) sprintf ('%s - %s', ... 
                                NAMES{pairs(i,:)}), (1 : size (pairs,1))', ...
                                'UniformOutput', false);
          NAMES = STATS.name;

      end

    end

    % Reorder fields in the STATS structure
    STATS.estimate = STATS.original;
    STATS = rmfield (STATS, 'original');
    switch (lower (METHOD))
      case 'wild'
        STATS.method = 'Wild bootstrap-t';
      case {'bayes','bayesian'}
        STATS.method = 'Bayesian bootstrap';
    end
    STATS = orderfields (STATS, {'method','name', 'estimate', 'CI_lower', ...
                                 'CI_upper', 'pval', 'fpr', 'n', 'prior'});

    % Print table of model coefficients and make figure of diagnostic plots
    switch (lower (DISPLAY))

      case {'on', true}

        % Print model formula
        fprintf('\nMODEL FORMULA (based on Wilkinson''s notation):\n\n%s\n', ...
                formula{end});

        % If applicable, print parameter estimates (a.k.a contrasts) for fixed
        % effects. Parameter estimates correspond to the contrasts we set.
        if (isempty (DIM))
          fprintf ('\nMODEL COEFFICIENTS\n\n');
          fprintf (cat (2, 'name                                   coeff', ...
                           '       CI_lower    CI_upper    p-val\n'));
          fprintf (cat (2, '--------------------------------------------', ...
                           '------------------------------------\n'));
        else
          switch (lower (POSTHOC))
            case 'none'
              fprintf ('\nMODEL ESTIMATED MARGINAL MEANS\n\n');
              fprintf (cat (2, 'name                                   ', ...
                               'mean        CI_lower    CI_upper        n\n'));
              fprintf (cat (2, '---------------------------------------', ...
                               '-----------------------------------------\n'));
            case {'pairwise', 'trt_vs_ctrl'}
              fprintf ('\nMODEL POSTHOC COMPARISONS\n\n');
              fprintf (cat (2, 'name                                   ', ...
                               'mean        CI_lower    CI_upper    p-adj\n'));
              fprintf (cat (2, '---------------------------------------', ...
                               '-----------------------------------------\n'));
          end
        end
        for j = 1:size (NAMES, 1)
          if ( (isempty (DIM)) || (ismember (lower (POSTHOC), ...
                                             {'pairwise', 'trt_vs_ctrl'})) )
            fprintf ('%-37s  %#-+10.4g  %#-+10.4g  %#-+10.4g', ...
                     NAMES{j}(1:min(end,37)), STATS.estimate(j), ...
                     STATS.CI_lower(j), STATS.CI_upper(j));
            if (isempty (STATS.pval))
              fprintf ('       \n');
            elseif (STATS.pval(j) <= 0.001)
              fprintf ('  <.001\n');
            elseif (STATS.pval(j) < 0.9995)
              fprintf ('   .%03u\n', round (STATS.pval(j) * 1e+03));
            elseif (isnan (STATS.pval(j)))
              fprintf ('    NaN\n');
            else
              fprintf ('  1.000\n');
            end
          else
            fprintf ('%-37s  %#-+10.4g  %#-+10.4g  %#-+10.4g  %5u\n', ...
                     NAMES{j}(1:min(end,37)), STATS.estimate(j), ...
                     STATS.CI_lower(j), STATS.CI_upper(j), STATS.n(j));
          end
        end
        fprintf('\n');

        % Make figure of diagnostic plots
        fhandle = figure (1);
        set (fhandle, 'Name', 'Diagnostic Plots: Model Residuals');
        h = diag (hat);                          % Leverage values
        mse = sse / dfe;                         % Mean squared error
        t = resid ./ (sqrt (mse * (1 - h)));     % Studentized residuals
        p = n - dfe;                             % Number of parameters
        fit = X * b;                             % Fitted values
        D = (1 / p) * t.^2 .* (h ./ (1 - h));    % Cook's distances
        [jnk, DI] = sort (D, 'descend');         % Sorted Cook's distances
        nk = 4;                                  % Number of most influential
                                                 % data points to label

        % Normal quantile-quantile plot
        subplot (2, 2, 1);
        x = ((1:n)' - .5) / n;
        [ts, I] = sort (t);
        stdnorminv = @(p) sqrt (2) * erfinv (2 * p - 1);
        q = stdnorminv (x);
        plot (q, ts, 'ok', 'markersize', 3);
        box off;
        grid on;
        xlabel ('Theoretical quantiles');
        ylabel ('Studentized Residuals');
        title ('Normal Q-Q Plot');
        arrayfun (@(i) text (q(I == DI(i)), t(DI(i)), ...
                             sprintf ('  %u', DI(i))), 1:min(nk,n))
        iqr = [0.25; 0.75]; 
        [ts, F] = bootcdf (t, true, 1);
        yl = interp1 (F, ts, iqr, 'linear', min (ts));
        xl = stdnorminv (iqr);
        slope = diff (yl) / diff (xl);
        int = yl(1) - slope * xl(1);
        ax1_xlim = get (gca, 'XLim');
        hold on; plot (ax1_xlim, slope * ax1_xlim + int, 'k-'); hold off;
        set (gca, 'Xlim', ax1_xlim);

        % Spread-Location Plot
        subplot (2, 2, 2);
        plot (fit, sqrt (abs (t)), 'ko', 'markersize', 3);
        box off;
        xlabel ('Fitted values');
        ylabel ('sqrt ( | Studentized Residuals | )');
        title ('Spread-Location Plot')
        ax2_xlim = get (gca, 'XLim');
        hold on; 
        plot (ax2_xlim, ones (1, 2) * sqrt (2), 'k:');
        plot (ax2_xlim, ones (1, 2) * sqrt (3), 'k-.'); 
        plot (ax2_xlim, ones (1, 2) * sqrt (4), 'k--');
        hold off;
        arrayfun (@(i) text (fit(DI(i)), sqrt (abs (t(DI(i)))), ...
                             sprintf ('  %u', DI(i))), [1:min(nk,n)]);
        xlim (ax2_xlim); 

        % Residual-Leverage plot
        subplot (2, 2, 3);
        plot (h, t, 'ko', 'markersize', 3);
        box off;
        xlabel ('Leverage')
        ylabel ('Studentized Residuals');
        title ('Residual-Leverage Plot')
        ax3_xlim = get (gca, 'XLim');
        ax3_ylim = get (gca, 'YLim');
        hold on; plot (ax3_xlim, zeros (1, 2), 'k-'); hold off;
        arrayfun (@(i) text (h(DI(i)), t(DI(i)), ...
                             sprintf ('  %u', DI(i))), [1:min(nk,n)]);
        set (gca, 'ygrid', 'on');
        xlim (ax3_xlim); ylim (ax3_ylim);

        % Cook's distance stem plot
        subplot (2, 2, 4);
        stem (D, 'ko', 'markersize', 3);
        box off;
        xlabel ('Obs. number')
        ylabel ('Cook''s distance')
        title ('Cook''s Distance Stem Plot')
        xlim ([0, n]);
        ax4_xlim = get (gca, 'XLim');
        ax4_ylim = get (gca, 'YLim');
        hold on; 
        plot (ax4_xlim, ones (1, 2) * 4 / dfe, 'k:');
        plot (ax4_xlim, ones (1, 2) * 0.5, 'k-.');
        plot (ax4_xlim, ones (1, 2), 'k--');
        hold off;
        arrayfun (@(i) text (DI(i), D(DI(i)), ...
                             sprintf ('  %u', DI(i))), [1:min(nk,n)]);
        xlim (ax4_xlim); ylim (ax4_ylim);

        set (findall ( gcf, '-property', 'FontSize'), 'FontSize', 7)

      case {'off', false}

        % do nothing

      otherwise

        error ('bootlm: wrong value for ''display'' parameter.')

  end

end

%--------------------------------------------------------------------------

function  [X, levels, nlevels, df, coeffnames, gid, CONTRASTS, ...
           center_continuous] = mDesignMatrix (GROUP, TERMS, ...
           CONTINUOUS, CONTRASTS, VARNAMES, n, Nm, Nx, Ng, cont_vec)

  % EVALUATE PREDICTOR LEVELS
  levels = cell (Nm, 1);
  gid = zeros (n, Nm);
  nlevels = zeros (Nm, 1);
  df = zeros (Nm + Nx, 1);
  termcols = ones (1 + Nm + Nx, 1);
  for j = 1:Nm
    if (any (j == CONTINUOUS))

      % CONTINUOUS PREDICTOR
      nlevels(j) = 1;
      termcols(j+1) = 1;
      df(j) = 1;
      if iscell (GROUP(:,j))
        gid(:,j) = cell2mat ([GROUP(:,j)]);
      else
        gid(:,j) = GROUP(:,j);
      end

    else

      % CATEGORICAL PREDICTOR
      levels{j} = unique_stable (GROUP(:,j));
      if isnumeric (levels{j})
        levels{j} = num2cell (levels{j});
      end
      nlevels(j) = numel (levels{j});
      for k = 1:nlevels(j)
        gid(ismember (GROUP(:,j),levels{j}{k}),j) = k;
      end
      termcols(j+1) = nlevels(j);
      df(j) = nlevels(j) - 1;

    end
  end

  % MAKE DESIGN MATRIX

  % MAIN EFFECTS
  X = cell (1, 1 + Nm + Nx);
  X{1} = ones (n, 1);
  coeffnames = cell (1, 1 + Nm + Nx);
  coeffnames{1} = '(Intercept)';
  vmeans = zeros (Nm, 1);
  center_continuous = cont_vec;
  for j = 1:Nm
    if (any (j == CONTINUOUS))

      % CONTINUOUS PREDICTOR
      if iscell (GROUP(:,j))
        X{1+j} = cell2mat (GROUP(:,j));
      else
        X{1+j} = GROUP(:,j);
      end
      if (strcmpi (CONTRASTS{j}, 'treatment'))
        % Don't center continuous variables if contrasts are 'treatment'
        center_continuous(j) = false;
        CONTRASTS{j} = [];
      else
        center_continuous(j) = true;
        vmeans(j) = mean ([X{1+j}]);
        X{1+j} = [X{1+j}] - vmeans(j);
      end
      % Create names of the coefficients relating to continuous main effects
      coeffnames{1+j} = VARNAMES{j};

    else

      % CATEGORICAL PREDICTOR
      if (isempty (CONTRASTS{j}))
        CONTRASTS{j} = contr_simple (nlevels(j));
      elseif (isnumeric (CONTRASTS{j}))
        % EVALUATE CUSTOM CONTRAST MATRIX
        % Check that the contrast matrix provided is the correct size
        if (~ all (size (CONTRASTS{j},1) == nlevels(j)))
          error (cat (2, 'bootlm: the number of rows in the contrast', ...
                         ' matrices should equal the number of', ...
                         ' predictor levels'))
        end
        if (~ all (size (CONTRASTS{j},2) == df(j)))
          error (cat (2, 'bootlm: the number of columns in each contrast', ...
                         ' matrix should equal the degrees of freedom', ...
                         ' (i.e. number of levels minus 1) for that predictor'))
        end
        if (~ all (any (CONTRASTS{j})))
          error (cat (2, 'bootlm: a contrast must be coded in each', ...
                         ' column of the contrast matrices'))
        end
      else
        switch (lower (CONTRASTS{j}))
          case {'simple','anova'}
            % SIMPLE EFFECT CODING (DEFAULT)
            % The first level is the reference level
            CONTRASTS{j} = contr_simple (nlevels(j));
          case 'poly'
            % POLYNOMIAL CONTRAST CODING
            CONTRASTS{j} = contr_poly (nlevels(j));
          case 'helmert'
            % HELMERT CONTRAST CODING
            CONTRASTS{j} = contr_helmert (nlevels(j));
          case 'effect'
            % DEVIATION EFFECT CONTRAST CODING
            CONTRASTS{j} = contr_sum (nlevels(j));
          case {'sdif','sdiff'}
            % SUCCESSIVE DEVIATIONS CONTRAST CODING
            CONTRASTS{j} = contr_sdif (nlevels(j));
          case 'treatment'
            % TREATMENT CONTRAST CODING
            % The first level is the reference level
            CONTRASTS{j} = contr_treatment (nlevels(j));
        end
      end
      C = CONTRASTS{j};
      func = @(x) x(gid(:,j));
      X{1+j} = cell2mat (cellfun (func, num2cell (C, 1), ...
                                  'UniformOutput', false));
      % Create names of the coefficients relating to continuous main effects
      coeffnames{1+j} = cell (df(j), 1);
      for v = 1:df(j)
        coeffnames{1+j}{v} = sprintf ('%s_%u', VARNAMES{j}, v);
      end

    end
  end

  % INTERACTION TERMS
  if (Nx > 0)
    row = TERMS((Ng > 1),:);
    for i = 1:Nx
      I = 1 + find (row(i,:));
      df(Nm+i) = prod (df(I-1));
      termcols(1+Nm+i) = prod (df(I-1) + 1);
      tmp = ones (n,1);
      for j = 1:numel(I);
        tmp = num2cell (tmp, 1);
        for k = 1:numel(tmp)
          tmp{k} = bsxfun (@times, tmp{k}, X{I(j)});
        end
        tmp = cell2mat (tmp);
      end
      X{1+Nm+i} = tmp;
      coeffnames{1+Nm+i} = cell (df(Nm+i),1);
      for v = 1:df(Nm+i)
        str = sprintf ('%s:', VARNAMES{I-1});
        coeffnames{1+Nm+i}{v} = strcat (str(1:end-1), '_', num2str (v));
      end
    end
  end

  % Remove any empty cells
  X = X(~ cellfun ('isempty', X));

end

%--------------------------------------------------------------------------

% BUILT IN CONTRAST CODING FUNCTIONS

function C = contr_simple (N)

  % Create contrast matrix (of doubles) using simple (ANOVA) contrast coding
  % These contrasts are centered (i.e. sum to 0)
  % Ideal for unordered predictors, with comparison to a reference level
  % The first predictor level is the reference level
  C =  cat (1, zeros (1, N - 1), eye (N - 1)) - (1 / N);

end

function C = contr_poly (N)

  % Create contrast matrix (of doubles) using polynomial contrast coding
  % for trend analysis of ordered categorical predictor levels
  % These contrasts are orthogonal and centered (i.e. sum to 0)
  % Ideal for ordered predictors
  [C, jnk] = qr (bsxfun (@power, (1:N)' - mean ((1:N)'), (0:(N - 1))));
  C(:,1) = [];
  s = ones (1, N - 1);
  s(1:2:N - 1) = s(1:2:N - 1) * -1;
  f = (sign(C(1,:)) ~= s);
  C(:,f) = C(:,f) * -1;

end

function C = contr_helmert (N)

  % Create contrast matrix (of doubles) using Helmert coding contrasts
  % These contrasts are orthogonal and centered (i.e. sum to 0)
  C = cat (1, tril (- ones (N - 1), -1) + diag ((N - 1):-1:1), ...
              -ones (1, N - 1)) ./ (N:-1:2);

end

function C = contr_sum (N)

  % Create contrast matrix (of doubles) using deviation effect coding
  % These contrasts are centered (i.e. sum to 0)
  C =  cat (1, - (ones (1, N - 1)), eye (N - 1));

end

function C = contr_sdif (N)

  % Create contrast matrix (of doubles) using successive differences coding
  % These contrasts are centered (i.e. sum to 0)
  C =  tril (ones (N, N - 1), -1) - ones (N, 1) / N * ((N - 1):-1:1);

end

function C = contr_treatment (N)

  % Create contrast matrix (of doubles) using treatment contrast coding
  % Ideal for unordered predictors, with comparison to a reference level
  % The first predictor level is the reference level
  C =  cat (1, zeros (1, N - 1), eye (N - 1));

end

%--------------------------------------------------------------------------

% FUNCTION TO FIT THE LINEAR MODEL

function [b, sse, resid, ucov, hat] = lmfit (X, Y, ISOCTAVE)

  % Get model coefficients by solving the linear equation. The number of free
  % parameters (i.e. intercept + coefficients) is equal to n - dfe (i.e. the
  % number of columns in X).
  b = X \ Y;                 % Equivalent to inv (X' * X) * (X' * y);

  % Get fitted values
  fit = X * b;

  % Get residuals from the fit
  resid = Y - fit;

  % Calculate the residual sums-of-squares
  sse = sum (resid.^2);

  % Calculate the unscaled covariance matrix (i.e. inv (X'*X )) and the Hat
  % matrix (i.e. X*(X'*X)^−1*X') by QR decomposition
  if (nargout > 3)
    [Q, R] = qr (X, 0);      % Economy-sized QR decomposition
    if ISOCTAVE
      ucov = chol2inv (R);
    else
      ucov = inv (R' * R);
    end
    hat = Q * Q';
  end

end

%--------------------------------------------------------------------------

% BUILT IN POSTHOC HYPOTHESIS TEST FUNCTIONS

function [L, pairs] = pairwise (L_EMM)

  % Get number of group members from the hypothesis matrix used 
  % to generate estimated marginal means
  Ng = size (unique (L_EMM', 'rows'), 1);

  % Create pairs matrix for pairwise comparisons
  gid = (1 : Ng)';  % Create numeric group ID
  A = ones (Ng, 1) * gid';
  B = tril (gid * ones(1, Ng),-1);
  pairs = [A(:), B(:)];
  ridx = (pairs(:, 2) == 0);
  pairs(ridx, :) = [];

  % Calculate hypothesis matrix for pairwise comparisons from the
  % estimated marginal means
  Np = size (pairs, 1);
  L_PWC = zeros (Np, Ng);
  for j = 1:Np
    L_PWC(j, pairs(j,:)) = [1,-1];
  end
  
  % Create hypothesis matrix to generate pairwise comparisons directly
  % from the regression coefficients. Note that the contrasts used to
  % fit the original model must sum to zero
  L = (L_PWC * L_EMM')';

end

%--------------------------------------------------------------------------

function [L, pairs] = trt_vs_ctrl (L_EMM, REF)

  if (nargin < 2)
    REF = 1;
  end

  % Get number of group members from the hypothesis matrix used 
  % to generate estimated marginal means
  Ng = size (unique (L_EMM','rows'), 1);
  if (REF > Ng)
    error ('trt_vs_ctrl: REF exceeds number of groups (i.e. rows in L_EMM)')
  end

  % Create pairs matrix for pairwise comparisons
  gid = (1 : Ng)';  % Create numeric group ID
  pairs = zeros (Ng - 1, 2);
  pairs(:, 1) = REF;
  pairs(:, 2) = gid(gid ~= REF);

  % Calculate hypothesis matrix for pairwise comparisons from the
  % estimated marginal means
  Np = size (pairs, 1);
  L_PWC = zeros (Np, Ng);
  for j = 1:Np
    L_PWC(j, pairs(j,:)) = [1,-1];
  end
  
  % Create hypothesis matrix to generate pairwise comparisons directly
  % from the regression coefficients. Note that the contrasts used to
  % fit the original model must sum to zero
  L = (L_PWC * L_EMM')';

end

%--------------------------------------------------------------------------

% FUNCTION TO CONTROL TYPE 1 ERROR ACROSS MULTIPLE POSTHOC COMPARISONS

function padj = holm (p)

  % Holm-Bonferroni procedure

  % Order raw p-values
  [ps, idx] = sort (p, 'ascend');
  k = numel (ps);

  % Implement Holm's step-down Bonferroni procedure
  padj = nan (k,1);
  padj(1) = k * ps(1);
  for j = 2:k
    padj(j) = max (padj(j - 1), (k - j + 1) * ps(j));
  end

  % Reorder the adjusted p-values to match the order of the original p-values
  [jnk, original_order] = sort (idx, 'ascend');
  padj = padj(original_order);

  % Truncate adjusted p-values to 1.0
  padj(padj>1) = 1;

end

%--------------------------------------------------------------------------

% FUNCTION TO FLATTEN A STRUCTURE ARRAY

function F = flatten_struct (S)

  fn = fieldnames (S);
  nm = numel (fn);
  F = struct;
  for i = 1:nm
    F.(fn{i}) = [S.(fn{i})]';
  end

end

%--------------------------------------------------------------------------

% FUNCTION THAT RETURNS UNIQUE VALUES IN THE ORDER THAT THEY FIRST APPEAR

function [U, IA, IC] = unique_stable (A, varargin)

  % Subfunction used for backwards compatibility

  % Error checking
  if any (ismember (varargin, {'first', 'last', 'sorted', 'stable'}))
    error ('unique_stable: the only option available is ''rows''')
  end
  if (iscell (A) && ismember ('rows', varargin))
    error ('unique_stable: ''rows'' option not supported for cell arrays')
  end

  % Flatten A to a column vector if 'rows' option is not specified
  if (~ ismember ('rows', varargin))
    A = A(:);
  end

  % Obtain sorted unique values
  [u, ia, ic] = unique (A, 'first', varargin{:});

  % Sort index of first occurence of unique values as they first appear
  IA = sort (ia);

  % Get unique values in the order of appearace (a.k.a. 'stable')
  U = A(IA,:);

  % Create vector of numeric identifiers for unique values in A
  n = numel (IA);
  if iscell (A)
    IC = sum (cell2mat (arrayfun (@(i) i * ismember (A, U(i,:)), ...
                        (1:n), 'UniformOutput', false)), 2);
  elseif isnumeric (A)
    IC = sum (cell2mat (arrayfun (@(i) i * (all (bsxfun (@eq, A, U(i,:)), ...
                        2)), (1:n), 'UniformOutput', false)), 2);
  end

end

%--------------------------------------------------------------------------

% FUNCTION TO PERFORM ANOVA

function AOVSTAT = bootanova (Y, X, DF, DFE, DEP, NBOOT, ALPHA, SEED, ISOCTAVE)

  % Bootstrap ANOVA (using sequential sums-of-squares, a.k.a. Type 1)

  % Compute observed statistics
  Nt = numel (DF) - 1;
  [jnk, SSE, RESID] = arrayfun (@(j) lmfit (X(:, 1 : sum (DF(1:j))), Y, ...
                                ISOCTAVE), (1:Nt + 1)', 'UniformOutput', false);
  SS = max (-diff (cell2mat (SSE)), 0);
  MS = SS ./ DF(2:end);
  MSE = SSE{end} / DFE;
  F = MS / MSE;

  % Obtain the F distribution under the null hypothesis by bootstrap of the
  % residuals from the full model. See ter Braak (1992) Permutation versus
  % bootstrap significance test in multiple regression and ANOVA. In Jockel
  % et al (Eds.) Bootstrapping and Related Techniques. Springer-Verlag, Berlin,
  % pg 79-86
  % See also the R function: https://rdrr.io/cran/lmboot/src/R/ANOVA.boot.R
  [jnk, jnk, BOOTSSE] = arrayfun (@(j) bootwild (RESID{end}, ...
                                X(:, 1 : sum (DF(1:j))), ...
                                DEP, NBOOT, ALPHA, SEED, [], ISOCTAVE), ...
                                (1:Nt + 1)', 'UniformOutput', false);
  BOOTSSE = cell2mat (BOOTSSE);
  BOOTSS = max (-diff (BOOTSSE), 0);
  BOOTMS = bsxfun (@rdivide, BOOTSS, DF(2:end));
  BOOTMSE = BOOTSSE(end,:) / DFE;
  BOOTF = bsxfun (@rdivide, BOOTMS, BOOTMSE);

  % Compute p-values
  res_lim = 1 / (NBOOT + 1);
  PVAL = nan (Nt, 1);
  for j = 1:Nt
    [x, jnk, P] = bootcdf (BOOTF(j,:), true, 1);
    if (F(j) < x(1))
      PVAL(j) = interp1 (x, P, F(j), 'linear', 1);
    else
      PVAL(j) = interp1 (x, P, F(j), 'linear', res_lim);
    end
  end

  % Prepare output
  AOVSTAT = struct ('MODEL', [], 'SS', SS, 'DF', DF(2:end), 'MS', MS, 'F', ...
                     F, 'PVAL', PVAL, 'SSE', SSE{end}, 'DFE', DFE, 'MSE', MSE);

end

%--------------------------------------------------------------------------

%!demo
%!
%! # Two-sample unpaired test on independent samples (equivalent to Welch's
%! # t-test). 
%!
%! score = [54 23 45 54 45 43 34 65 77 46 65]';
%! gender = {'male' 'male' 'male' 'male' 'male' 'female' 'female' 'female' ...
%!           'female' 'female' 'female'}';
%!
%! % 95% confidence intervals and p-values for the difference in mean score
%! % between males and females (computed by wild bootstrap)
%! STATS = bootlm (score, gender, 'display', 'on', 'varnames', 'gender', ...
%!                 'dim', 1, 'posthoc','trt_vs_ctrl');
%!
%! % 95% credible intervals for the estimated marginal means of the scores by
%! % males and females (computed by Bayesian bootstrap)
%! STATS = bootlm (score, gender, 'display', 'on', 'varnames', 'gender', ...
%!                 'dim', 1, 'method', 'bayesian', 'prior', 'auto');


%!demo
%!
%! # Two-sample paired test on dependent or matched samples equivalent to a
%! # paired t-test.
%!
%! score = [4.5 5.6; 3.7 6.4; 5.3 6.4; 5.4 6.0; 3.9 5.7]';
%! treatment = {'before' 'after'; 'before' 'after'; 'before' 'after';
%!              'before' 'after'; 'before' 'after'}';
%! subject = {'GS' 'GS'; 'JM' 'JM'; 'HM' 'HM'; 'JW' 'JW'; 'PS' 'PS'}';
%!
%! % 95% confidence intervals and p-values for the difference in mean score
%! % before and after treatment (computed by wild bootstrap)
%! STATS = bootlm (score(:), {treatment(:), subject(:)}, ...
%!                            'model', 'linear', 'display', 'on', ...
%!                            'varnames', {'treatment', 'subject'}, ...
%!                            'dim', 1, 'posthoc','trt_vs_ctrl');
%!
%! % 95% credible intervals for the estimated marginal means of the scores
%! % before and after treatment (computed by Bayesian bootstrap)
%! STATS = bootlm (score(:), {treatment(:), subject(:)}, ...
%!                            'model', 'linear', 'display', 'on', ...
%!                            'varnames', {'treatment', 'subject'}, ...
%!                            'dim', 1, 'method','bayesian', 'prior', 'auto');

%!demo
%!
%! # One-way design. The data is from a study on the strength of structural
%! # beams, in Hogg and Ledolter (1987) Engineering Statistics. NY: MacMillan
%!
%! strength = [82 86 79 83 84 85 86 87 74 82 ...
%!            78 75 76 77 79 79 77 78 82 79]';
%! alloy = {'st','st','st','st','st','st','st','st', ...
%!          'al1','al1','al1','al1','al1','al1', ...
%!          'al2','al2','al2','al2','al2','al2'}';
%!
%! % 95% confidence intervals and p-values for the differences in mean strength
%! % of three alloys (computed by wild bootstrap)
%! STATS = bootlm (strength, alloy, 'display', 'on', 'varnames', 'alloy', ...
%!                 'dim', 1, 'posthoc','pairwise');
%!
%! % 95% credible intervals for the estimated marginal means of the strengths
%! % of each of the alloys (computed by Bayesian bootstrap)
%! STATS = bootlm (strength, alloy, 'display', 'on', 'varnames', 'alloy', ...
%!                 'dim', 1, 'method','bayesian', 'prior', 'auto');

%!demo
%!
%! # One-way repeated measures design. The data is from a study on the number of
%! # words recalled by 10 subjects for three time condtions, in Loftus & Masson
%! # (1994) Psychon Bull Rev. 1(4):476-490, Table 2.
%!
%! words = [10 13 13; 6 8 8; 11 14 14; 22 23 25; 16 18 20; ...
%!          15 17 17; 1 1 4; 12 15 17;  9 12 12;  8 9 12];
%! seconds = [1 2 5; 1 2 5; 1 2 5; 1 2 5; 1 2 5; ...
%!            1 2 5; 1 2 5; 1 2 5; 1 2 5; 1 2 5;];
%! subject = [ 1  1  1;  2  2  2;  3  3  3;  4  4  4;  5  5  5; ...
%!             6  6  6;  7  7  7;  8  8  8;  9  9  9; 10 10 10];
%!
%! % 95% confidence intervals and p-values for the differences in mean number of
%! % words recalled for the different times (using wild bootstrap).
%! STATS = bootlm (words(:), {seconds(:), subject(:)}, ...
%!                            'model', 'linear', 'display', 'on', ...
%!                            'varnames', {'seconds', 'subject'}, ...
%!                            'dim', 1, 'posthoc', 'pairwise');
%!
%! % 95% credible intervals for the estimated marginal means of the number of
%! % words recalled for each time (computed using Bayesian bootstrap).
%! STATS = bootlm (words(:), {seconds(:), subject(:)}, ...
%!                            'model', 'linear', 'display', 'on', ...
%!                            'varnames', {'seconds', 'subject'}, ...
%!                            'dim', 1, 'method', 'bayesian', 'prior', 'auto');

%!demo
%!
%! # Balanced two-way design. The data is yield of cups of popped popcorn from
%! # different popcorn brands and popper types, in Hogg and Ledolter (1987)
%! # Engineering Statistics. NY: MacMillan
%!
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! brands = {'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'};
%! popper = {'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; ...
%!           'air', 'air', 'air'; 'air', 'air', 'air'; 'air', 'air', 'air'};
%!
%! % Check regression coefficients corresponding to brand x popper interaction
%! % using 'anova' contrast coding
%! STATS = bootlm (popcorn(:), {brands(:), popper(:)}, ...
%!                            'display', 'on', 'model', 'full', ...
%!                            'varnames', {'brands', 'popper'});
%!
%! % 95% confidence intervals and p-values for the differences in mean yield of
%! % different popcorn brands (computed by wild bootstrap).
%! STATS = bootlm (popcorn(:), {brands(:), popper(:)}, ...
%!                            'display', 'on', 'model', 'full', ...
%!                            'varnames', {'brands', 'popper'}, ...
%!                            'dim', 1, 'posthoc', 'pairwise');
%!
%! % 95% credible intervals for the estimated marginal means of the yield for
%! % each popcorn brand (computed by Bayesian bootstrap).
%! STATS = bootlm (popcorn(:), {brands(:), popper(:)}, ...
%!                            'display', 'on', 'model', 'full', ...
%!                            'varnames', {'brands', 'popper'}, ...
%!                            'dim', 1, 'method', 'bayesian', 'prior', 'auto');
%!
%! % 95% confidence intervals and p-values for the differences in mean yield
%! % for different popper types (computed by wild bootstrap).
%! STATS = bootlm (popcorn(:), {brands(:), popper(:)}, ...
%!                            'display', 'on', 'model', 'full', ...
%!                            'varnames', {'brands', 'popper'}, ...
%!                            'dim', 2, 'posthoc', 'pairwise');
%!
%! % 95% credible intervals for the estimated marginal means of the yield for
%! % each popper type (computed by Bayesian bootstrap).
%! STATS = bootlm (popcorn(:), {brands(:), popper(:)}, ...
%!                            'display', 'on', 'model', 'full', ...
%!                            'varnames', {'brands', 'popper'}, ...
%!                            'dim', 2, 'method', 'bayesian', 'prior', 'auto');

%!demo
%!
%! # Unbalanced two-way design (2x2). The data is from a study on the effects
%! # of gender and a college degree on starting salaries of company employees,
%! # in Maxwell, Delaney and Kelly (2018): Chapter 7, Table 15
%!
%! salary = [24 26 25 24 27 24 27 23 15 17 20 16, ...
%!           25 29 27 19 18 21 20 21 22 19]';
%! gender = {'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f'...
%!           'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm'}';
%! degree = [1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 0 0 0 0 0 0 0]';
%!
%! % ANOVA (including the main effect of gender averaged over levels of degree)
%! [STATS, BOOTSTAT, AOVSTAT] = bootlm (salary, {degree, gender}, 'model', ...
%!                             'full', 'display', 'off', 'varnames', ...
%!                             {'degree', 'gender'});
%!
%! fprintf ('ANOVA SUMMARY with gender averaged over levels of degree\n')
%! for i = 1:numel(AOVSTAT.F)
%!   fprintf ('F(%u,%u) = %.2f, p = %.3g for the model: %s\n', ...
%!            AOVSTAT.DF(i), AOVSTAT.DFE, AOVSTAT.F(i), ...
%!            AOVSTAT.PVAL(i), AOVSTAT.MODEL{i});
%! end
%!
%! % ANOVA (including the main effect of degree averaged over levels of gender)
%! [STATS, BOOTSTAT, AOVSTAT] = bootlm (salary, {gender, degree}, 'model', ...
%!                             'full', 'display', 'off', 'varnames', ...
%!                             {'gender', 'degree'});
%!
%! fprintf ('\nANOVA SUMMARY with degree averaged over levels of gender\n')
%! for i = 1:numel(AOVSTAT.F)
%!   fprintf ('F(%u,%u) = %.2f, p = %.3g for the model: %s\n', ...
%!            AOVSTAT.DF(i), AOVSTAT.DFE, AOVSTAT.F(i), ...
%!            AOVSTAT.PVAL(i), AOVSTAT.MODEL{i});
%! end
%!
%! % Check regression coefficient corresponding to gender x degree interaction
%! % using 'anova' contrast coding
%! STATS = bootlm (salary, {gender, degree}, 'model', 'full', ...
%!                             'display', 'on', 'varnames', ...
%!                             {'gender', 'degree'});
%!
%! % 95% confidence intervals and p-values for the differences in mean salary
%! % between males and females (computed by wild bootstrap).
%! STATS = bootlm (salary, {gender, degree}, 'model', 'full', ...
%!                            'display', 'on', 'varnames', ...
%!                            {'gender', 'degree'}, 'dim', 1, ...
%!                            'posthoc', 'trt_vs_ctrl');
%!
%! % 95% credible intervals for the estimated marginal means for salaries of
%! % females and males (computed by Bayesian bootstrap).
%! STATS = bootlm (salary, {gender, degree}, 'model', 'full', ...
%!                            'display', 'on', 'varnames', ...
%!                            {'gender', 'degree'}, 'dim', 1, ...
%!                            'method', 'bayesian', 'prior', 'auto');
%!
%! % 95% confidence intervals and p-values for the differences in mean salary
%! % between employees with or without a degree (computed by wild bootstrap).
%! STATS = bootlm (salary, {gender, degree}, 'model', 'full', ...
%!                            'display', 'on', 'varnames', ...
%!                            {'gender', 'degree'}, 'dim', 2, ...
%!                            'posthoc', 'trt_vs_ctrl');
%!
%! % 95% credible intervals for the estimated marginal means for salaries of
%! % employees with or without a degree (computed by Bayesian bootstrap).
%! STATS = bootlm (salary, {gender, degree}, 'model', 'full', ...
%!                            'display', 'on', 'varnames', ...
%!                            {'gender', 'degree'}, 'dim', 2, ...
%!                            'method', 'bayesian','prior', 'auto');

%!demo
%!
%! # Unbalanced three-way design (3x2x2). The data is from a study of the
%! # effects of three different drugs, biofeedback and diet on patient blood
%! # pressure, adapted* from Maxwell, Delaney and Kelly (2018): Ch 8, Table 12
%!
%! drug = {'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' ...
%!         'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X';
%!         'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' ...
%!         'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y';
%!         'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' ...
%!         'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z'};
%! feedback = [1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
%!             1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
%!             1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0];
%! diet = [0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1;
%!         0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1;
%!         0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1];
%! BP = [170 175 165 180 160 158 161 173 157 152 181 190 ...
%!       173 194 197 190 176 198 164 190 169 164 176 175;
%!       186 194 201 215 219 209 164 166 159 182 187 174 ...
%!       189 194 217 206 199 195 171 173 196 199 180 203;
%!       180 187 199 170 204 194 162 184 183 156 180 173 ...
%!       202 228 190 206 224 204 205 199 170 160 179 179];
%!
%! % Perform 3-way ANOVA (this design is balanced so order of predictors does 
%! % not make any difference)
%! [STATS, BOOTSTAT, AOVSTAT] = bootlm (BP(:), {drug(:), feedback(:), ...
%!                                    diet(:)}, 'seed', 1, ...
%!                                    'model', 'full', 'display', 'off', ...
%!                                    'varnames', {'drug', 'feedback', 'diet'});
%!
%! fprintf ('ANOVA SUMMARY\n')
%! for i = 1:numel(AOVSTAT.F)
%!   fprintf ('F(%u,%u) = %.2f, p = %.3g for the model: %s\n', ...
%!            AOVSTAT.DF(i), AOVSTAT.DFE, AOVSTAT.F(i), ...
%!            AOVSTAT.PVAL(i), AOVSTAT.MODEL{i});
%! end
%!
%! % Check regression coefficient corresponding to drug x feedback x diet
%! % interaction using 'anova' contrast coding
%! STATS = bootlm (BP(:), {drug(:), feedback(:), diet(:)}, ...
%!                                    'model', 'full', ...
%!                                    'display', 'on', ...
%!                                    'varnames', {'drug', 'feedback', 'diet'});
%!
%! % 95% confidence intervals and p-values for the differences in mean salary
%! % between males and females (computed by wild bootstrap).
%! STATS = bootlm (BP(:), {drug(:), feedback(:), diet(:)}, 'model', 'full', ...
%!                                    'display', 'on', 'dim', [1,2,3], ...
%!                                    'posthoc', 'trt_vs_ctrl', ...
%!                                    'varnames', {'drug', 'feedback', 'diet'});
%!
%! % 95% credible intervals for the estimated marginal means of salaries of
%! % females and males (computed by Bayesian bootstrap).
%! STATS = bootlm (BP(:), {drug(:), feedback(:), diet(:)}, 'model', 'full', ...
%!                                    'display', 'on', 'dim', [1,2,3], ...
%!                                    'method', 'bayesian', 'prior', 'auto', ...
%!                                    'varnames', {'drug', 'feedback', 'diet'});

%!demo
%!
%! # One-way design with continuous covariate. The data is from a study of the
%! # additive effects of species and temperature on chirpy pulses of crickets,
%! # from Stitch, The Worst Stats Text eveR
%!
%! pulse = [67.9 65.1 77.3 78.7 79.4 80.4 85.8 86.6 87.5 89.1 ...
%!          98.6 100.8 99.3 101.7 44.3 47.2 47.6 49.6 50.3 51.8 ...
%!          60 58.5 58.9 60.7 69.8 70.9 76.2 76.1 77 77.7 84.7]';
%! temp = [20.8 20.8 24 24 24 24 26.2 26.2 26.2 26.2 28.4 ...
%!         29 30.4 30.4 17.2 18.3 18.3 18.3 18.9 18.9 20.4 ...
%!         21 21 22.1 23.5 24.2 25.9 26.5 26.5 26.5 28.6]';
%! species = {'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' ...
%!            'ex' 'ex' 'ex' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' ...
%!            'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv'};
%!
%! % Perform ANCOVA 
%! [STATS, BOOTSTAT, AOVSTAT] = bootlm (pulse, {temp, species}, 'model', ...
%!                           'linear', 'continuous', 1, 'display', 'off', ...
%!                           'varnames', {'temp', 'species'});
%!
%! fprintf ('ANCOVA SUMMARY\n')
%! for i = 1:numel(AOVSTAT.F)
%!   fprintf ('F(%u,%u) = %.2f, p = %.3g for the model: %s\n', ...
%!            AOVSTAT.DF(i), AOVSTAT.DFE, AOVSTAT.F(i), ...
%!            AOVSTAT.PVAL(i), AOVSTAT.MODEL{i});
%! end
%!
%! % Estimate regression coefficients using 'anova' contrast coding 
%! STATS = bootlm (pulse, {temp, species}, 'model', 'linear', ...
%!                           'continuous', 1, 'display', 'on', ...
%!                           'varnames', {'temp', 'species'});
%!
%! % 95% confidence intervals and p-values for the differences in the mean of
%! % chirpy pulses of ex ad niv species (computed by wild bootstrap).
%! STATS = bootlm (pulse, {temp, species}, 'model', 'linear', ...
%!                           'continuous', 1, 'display', 'on', ...
%!                           'varnames', {'temp', 'species'}, 'dim', 2, ...
%!                           'posthoc', 'trt_vs_ctrl');
%!
%! % 95% credible intervals for the estimated marginal means of chirpy pulses
%! % of ex and niv species (computed by Bayesian bootstrap).
%! STATS = bootlm (pulse, {temp, species}, 'model', 'linear', ...
%!                           'continuous', 1, 'display', 'on', ...
%!                           'varnames', {'temp', 'species'}, 'dim', 2, ...
%!                           'method', 'bayesian', 'prior', 'auto');

%!demo
%!
%! # Factorial design with continuous covariate. The data is from a study of the
%! # effects of treatment and exercise on stress reduction score after adjusting
%! # for age. Data from R datarium package).
%!
%! score = [95.6 82.2 97.2 96.4 81.4 83.6 89.4 83.8 83.3 85.7 ...
%!          97.2 78.2 78.9 91.8 86.9 84.1 88.6 89.8 87.3 85.4 ...
%!          81.8 65.8 68.1 70.0 69.9 75.1 72.3 70.9 71.5 72.5 ...
%!          84.9 96.1 94.6 82.5 90.7 87.0 86.8 93.3 87.6 92.4 ...
%!          100. 80.5 92.9 84.0 88.4 91.1 85.7 91.3 92.3 87.9 ...
%!          91.7 88.6 75.8 75.7 75.3 82.4 80.1 86.0 81.8 82.5]';
%! treatment = {'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' ...
%!              'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' ...
%!              'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' ...
%!              'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  ...
%!              'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  ...
%!              'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'}';
%! exercise = {'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  ...
%!             'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' ...
%!             'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  ...
%!             'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  ...
%!             'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' ...
%!             'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'}';
%! age = [59 65 70 66 61 65 57 61 58 55 62 61 60 59 55 57 60 63 62 57 ...
%!        58 56 57 59 59 60 55 53 55 58 68 62 61 54 59 63 60 67 60 67 ...
%!        75 54 57 62 65 60 58 61 65 57 56 58 58 58 52 53 60 62 61 61]';
%!
%! % ANOVA/ANCOVA statistics
%! [STATS, BOOTSTAT, AOVSTAT] = bootlm (score, {age, exercise, treatment}, ...
%!                            'model', [1 0 0; 0 1 0; 0 0 1; 0 1 1], ...
%!                            'continuous', 1, 'display', 'off', ...
%!                            'varnames', {'age', 'exercise', 'treatment'});
%!
%! fprintf ('ANOVA / ANCOVA SUMMARY\n')
%! for i = 1:numel(AOVSTAT.F)
%!   fprintf ('F(%u,%u) = %.2f, p = %.3g for the model: %s\n', ...
%!            AOVSTAT.DF(i), AOVSTAT.DFE, AOVSTAT.F(i), ...
%!            AOVSTAT.PVAL(i), AOVSTAT.MODEL{i});
%! end
%!
%! % Estimate regression coefficients using 'anova' contrast coding 
%! STATS = bootlm (score, {age, exercise, treatment}, ...
%!                            'model', [1 0 0; 0 1 0; 0 0 1; 0 1 1], ...
%!                            'continuous', 1, 'display', 'on', ...
%!                            'varnames', {'age', 'exercise', 'treatment'});
%!
%! % 95% confidence intervals and p-values for the differences in mean score
%! % across different treatments and amounts of exercise after adjusting for
%  % age (computed by wild bootstrap).
%! STATS = bootlm (score, {age, exercise, treatment}, ...
%!                            'model', [1 0 0; 0 1 0; 0 0 1; 0 1 1], ...
%!                            'continuous', 1, 'display', 'on', ...
%!                            'varnames', {'age', 'exercise', 'treatment'}, ...
%!                            'dim', [2, 3], 'posthoc', 'trt_vs_ctrl');
%!
%! % 95% credible intervals for the estimated marginal means of scores across
%! % different treatments and amounts of exercise after adjusting for age
%! % (computed by Bayesian bootstrap).
%! STATS = bootlm (score, {age, exercise, treatment}, 'dim', [2, 3], ...
%!                            'model', [1 0 0; 0 1 0; 0 0 1; 0 1 1], ...
%!                            'continuous', 1, 'display', 'on', ...
%!                            'varnames', {'age', 'exercise', 'treatment'}, ...
%!                            'method', 'bayesian', 'prior', 'auto');

%!demo
%!
%! # Unbalanced one-way design with custom, orthogonal contrasts. The statistics
%! # relating to the contrasts are shown in the table of model parameters, and
%! # can be retrieved from the STATS.coeffs output.
%!
%! dv =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! g = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 ...
%!      4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%! C = [ 0.4001601  0.3333333  0.5  0.0
%!       0.4001601  0.3333333 -0.5  0.0
%!       0.4001601 -0.6666667  0.0  0.0
%!      -0.6002401  0.0000000  0.0  0.5
%!      -0.6002401  0.0000000  0.0 -0.5];
%!
%! % 95% confidence intervals and p-values for linear contrasts 
%! STATS = bootlm (dv, g, 'contrasts', C, 'varnames', 'score', ...
%!                          'alpha', 0.05, 'display', true);
%!
%! % 95% credible intervals for estimated marginal means 
%! STATS = bootlm (dv, g, 'contrasts', C, 'varnames', 'score', ...
%!                          'alpha', 0.05, 'display', true, 'dim', 1, ...
%!                          'method', 'Bayesian', 'prior', 'auto');


%!test
%!
%! # Two-sample unpaired test on independent samples (equivalent to Welch's
%! # t-test).
%!
%! score = [54 23 45 54 45 43 34 65 77 46 65]';
%! gender = {'male' 'male' 'male' 'male' 'male' 'female' 'female' 'female' ...
%!           'female' 'female' 'female'}';
%!
%! stats = bootlm (score, gender, 'display', 'off', 'varnames', 'gender', ...
%!                                'seed', 1);
%!
%! assert (stats.pval(2), 0.2434934955512797, 1e-09);
%! assert (stats.fpr(2), 0.4832095599189747, 1e-09);
%! # ttest2 (with 'vartype' = 'unequal') gives a p-value of 0.2501;

%!test
%!
%! # Two-sample paired test on dependent or matched samples equivalent to a
%! # paired t-test.
%!
%! score = [4.5 5.6; 3.7 6.4; 5.3 6.4; 5.4 6.0; 3.9 5.7]';
%! treatment = {'before' 'after'; 'before' 'after'; 'before' 'after';
%!              'before' 'after'; 'before' 'after'}';
%! subject = {'GS' 'GS'; 'JM' 'JM'; 'HM' 'HM'; 'JW' 'JW'; 'PS' 'PS'}';
%!
%! stats = bootlm (score(:), {treatment(:), subject(:)}, 'seed', 1, ...
%!                            'model', 'linear', 'display', 'off', ...
%!                            'varnames', {'treatment', 'subject'});
%!
%! assert (stats.pval(1), 0.0007121854921651461, 1e-09);
%! assert (stats.pval(2), 0.002663469844077049, 1e-09);
%! assert (stats.pval(3), 0.9999999999999917, 1e-09);
%! assert (stats.pval(4), 0.06635496003290851, 1e-09);
%! assert (stats.pval(5), 0.4382333666561282, 1e-09);
%! assert (stats.pval(6), 0.3639361232818474, 1e-09);

%!test
%!
%! # One-way design. The data is from a study on the strength of structural
%! # beams, in Hogg and Ledolter (1987) Engineering Statistics. NY: MacMillan
%!
%! strength = [82 86 79 83 84 85 86 87 74 82 ...
%!            78 75 76 77 79 79 77 78 82 79]';
%! alloy = {'st','st','st','st','st','st','st','st', ...
%!          'al1','al1','al1','al1','al1','al1', ...
%!          'al2','al2','al2','al2','al2','al2'}';
%!
%! stats = bootlm (strength, alloy, 'display', 'off', 'varnames', 'alloy', ...
%!                                  'seed', 1);
%!
%! assert (stats.CI_lower(2), -10.17909151307657, 1e-09);
%! assert (stats.CI_upper(2), -3.820908486923432, 1e-09);
%! assert (stats.CI_lower(3), -7.462255988161777, 1e-09);
%! assert (stats.CI_upper(3), -2.537744011838216, 1e-09);

%!test
%!
%! # One-way repeated measures design. The data is from a study on the number of
%! # words recalled by 10 subjects for three time condtions, in Loftus & Masson
%! # (1994) Psychon Bull Rev. 1(4):476-490, Table 2.
%!
%! words = [10 13 13; 6 8 8; 11 14 14; 22 23 25; 16 18 20; ...
%!          15 17 17; 1 1 4; 12 15 17;  9 12 12;  8 9 12];
%! seconds = [1 2 5; 1 2 5; 1 2 5; 1 2 5; 1 2 5; ...
%!            1 2 5; 1 2 5; 1 2 5; 1 2 5; 1 2 5;];
%! subject = [ 1  1  1;  2  2  2;  3  3  3;  4  4  4;  5  5  5; ...
%!             6  6  6;  7  7  7;  8  8  8;  9  9  9; 10 10 10];
%!
%! stats = bootlm (words(:), {seconds(:), subject(:)}, 'seed', 1, ...
%!                            'model', 'linear', 'display', 'off', ...
%!                            'varnames', {'seconds', 'subject'});
%!
%! assert (stats.CI_lower(2), 1.266092224054235, 1e-09);
%! assert (stats.CI_upper(2), 2.733907775945761, 1e-09);
%! assert (stats.CI_lower(3), 2.554265809089302, 1e-09);
%! assert (stats.CI_upper(3), 3.845734190910699, 1e-09);

%!test
%!
%! # Balanced two-way design with interaction. The data is from a study of
%! # popcorn brands and popper types, in Hogg and Ledolter (1987) Engineering
%! # Statistics. New York: MacMillan
%!
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! brands = {'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'};
%! popper = {'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; ...
%!           'air', 'air', 'air'; 'air', 'air', 'air'; 'air', 'air', 'air'};
%!
%! stats = bootlm (popcorn(:), {brands(:), popper(:)}, 'seed', 1, ...
%!                            'display', 'off', 'model', 'full', ...
%!                            'varnames', {'brands', 'popper'});
%!
%! assert (stats.pval(2), 0.0001, 1e-09);
%! assert (stats.pval(3), 0.0001, 1e-09);
%! assert (stats.pval(4), 0.0001338515324481833, 1e-09);
%! assert (stats.pval(5), 0.3403340334033442, 1e-09);
%! assert (stats.pval(6), 0.7317882724305511, 1e-09);
%! assert (stats.fpr(2), 0.00249737757706675, 1e-09);
%! assert (stats.fpr(3), 0.00249737757706675, 1e-09);
%! assert (stats.fpr(4), 0.003234567489304454, 1e-09);
%! assert (stats.fpr(5), 0.4992799823055505, 1e-09);
%! assert (stats.fpr(6), 0.5, 1e-09);


%!test
%!
%! # Unbalanced two-way design (2x2). The data is from a study on the effects
%! # of gender and having a college degree on salaries of company employees,
%! # in Maxwell, Delaney and Kelly (2018): Chapter 7, Table 15
%!
%! salary = [24 26 25 24 27 24 27 23 15 17 20 16, ...
%!           25 29 27 19 18 21 20 21 22 19]';
%! gender = {'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f'...
%!           'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm'}';
%! degree = [1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 0 0 0 0 0 0 0]';
%!
%! [stats, bootstat, aovstats] = bootlm (salary, {gender, degree}, ...
%!                            'model', 'full', 'display', 'off', 'varnames', ...
%!                            {'gender', 'degree'}, 'seed', 1);
%!
%! assert (aovstats.PVAL(1), 0.7523035992551597, 1e-09);   % Normal ANOVA: 0.747 
%! assert (aovstats.PVAL(2), 0.0001, 1e-09);               % Normal ANOVA: <.001 
%! assert (aovstats.PVAL(3), 0.5666177238662272, 1e-09);   % Normal ANOVA: 0.524
%! assert (stats.pval(2), 0.01257386294526463, 1e-09);
%! assert (stats.pval(3), 0.0001, 1e-09);
%! assert (stats.pval(4), 0.5820694859231055, 1e-09);
%! assert (stats.fpr(2), 0.1301119743071732, 1e-09);
%! assert (stats.fpr(3), 0.00249737757706675, 1e-09);
%! assert (stats.fpr(4), 0.5, 1e-09);
%!
%! [stats, bootstat, aovstats] = bootlm (salary, {degree, gender}, ...
%!                            'model', 'full', 'display', 'off', 'varnames', ...
%!                            {'degree', 'gender'}, 'seed', 1);
%!
%! assert (aovstats.PVAL(1), 0.0001, 1e-09);               % Normal ANOVA: <.001 
%! assert (aovstats.PVAL(2), 0.004950446391560281, 1e-09); % Normal ANOVA: 0.004
%! assert (aovstats.PVAL(3), 0.566617723866227, 1e-09);    % Normal ANOVA: 0.524
%! assert (stats.pval(2), 0.0001, 1e-09);
%! assert (stats.pval(3), 0.01257386294526464, 1e-09);
%! assert (stats.pval(4), 0.5820694859231067, 1e-09);
%! assert (stats.fpr(2), 0.00249737757706675, 1e-09);
%! assert (stats.fpr(3), 0.1301119743071733, 1e-09);
%! assert (stats.fpr(4), 0.5, 1e-09);

%!test
%!
%! # Unbalanced two-way design (3x2). The data is from a study of the effect of
%! # adding sugar and/or milk on the tendency of coffee to make people babble,
%! # in from Navarro (2019): 16.10
%!
%! sugar = {'real' 'fake' 'fake' 'real' 'real' 'real' 'none' 'none' 'none' ...
%!          'fake' 'fake' 'fake' 'real' 'real' 'real' 'none' 'none' 'fake'}';
%! milk = {'yes' 'no' 'no' 'yes' 'yes' 'no' 'yes' 'yes' 'yes' ...
%!         'no' 'no' 'yes' 'no' 'no' 'no' 'no' 'no' 'yes'}';
%! babble = [4.6 4.4 3.9 5.6 5.1 5.5 3.9 3.5 3.7...
%!           5.6 4.7 5.9 6.0 5.4 6.6 5.8 5.3 5.7]';
%!
%! stats = bootlm (babble, {sugar, milk}, 'model', 'full', 'display', 'off', ...
%!                                'seed', 1, 'varnames', {'sugar', 'milk'});
%!
%! assert (stats.pval(5), 0.00433268463709287, 1e-09);
%! assert (stats.pval(6), 0.05620134119970051, 1e-09);
%! assert (stats.fpr(5), 0.06022795764518322, 1e-09);
%! assert (stats.fpr(6), 0.305458916611921, 1e-09);

%!test
%!
%! # Balanced three-way design (3x2x2). The data is from a study of the
%! # effects of three different drugs, biofeedback and diet on patient blood
%! # pressure, adapted* from Maxwell, Delaney and Kelly (2018): Ch 8, Table 12
%!
%! drug = {'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' ...
%!         'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X';
%!         'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' ...
%!         'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y';
%!         'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' ...
%!         'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z'};
%! feedback = [1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
%!             1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
%!             1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0];
%! diet = [0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1;
%!         0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1;
%!         0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1];
%! BP = [170 175 165 180 160 158 161 173 157 152 181 190 ...
%!       173 194 197 190 176 198 164 190 169 164 176 175;
%!       186 194 201 215 219 209 164 166 159 182 187 174 ...
%!       189 194 217 206 199 195 171 173 196 199 180 203;
%!       180 187 199 170 204 194 162 184 183 156 180 173 ...
%!       202 228 190 206 224 204 205 199 170 160 179 179];
%!
%! [stats, bootstat, aovstat] = bootlm (BP(:), {drug(:), feedback(:), ...
%!                                    diet(:)}, 'seed', 1, ...
%!                                    'model', 'full', 'display', 'off', ...
%!                                    'varnames', {'drug', 'feedback', 'diet'});
%!
%! assert (aovstat.PVAL(1), 0.0001785474921424416, 1e-09);
%! assert (aovstat.PVAL(2), 0.0005607720210921765, 1e-09);
%! assert (aovstat.PVAL(3), 0.0001, 1e-09);
%! assert (aovstat.PVAL(4), 0.4343155166545469, 1e-09);
%! assert (aovstat.PVAL(5), 0.06277877943312708, 1e-09);
%! assert (aovstat.PVAL(6), 0.6484269049223992, 1e-09);
%! assert (aovstat.PVAL(7), 0.03878235882690048, 1e-09);
%!
%! stats = bootlm (BP(:), {drug(:), feedback(:), diet(:)}, ...
%!                                    'seed', 1, ...
%!                                    'model', 'full', ...
%!                                    'display', 'off', ...
%!                                    'varnames', {'drug', 'feedback', 'diet'});
%!
%! assert (stats.pval(11), 0.01390969497190977, 1e-09);
%! assert (stats.pval(12), 0.7334357687702242, 1e-09);
%! assert (stats.fpr(11), 0.1391526669998761, 1e-09);
%! assert (stats.fpr(12), 0.5, 1e-09);

%!test
%!
%! # One-way design with continuous covariate. The data is from a study of the
%! # additive effects of species and temperature on chirpy pulses of crickets,
%! # from Stitch, The Worst Stats Text eveR
%!
%! pulse = [67.9 65.1 77.3 78.7 79.4 80.4 85.8 86.6 87.5 89.1 ...
%!          98.6 100.8 99.3 101.7 44.3 47.2 47.6 49.6 50.3 51.8 ...
%!          60 58.5 58.9 60.7 69.8 70.9 76.2 76.1 77 77.7 84.7]';
%! temp = [20.8 20.8 24 24 24 24 26.2 26.2 26.2 26.2 28.4 ...
%!         29 30.4 30.4 17.2 18.3 18.3 18.3 18.9 18.9 20.4 ...
%!         21 21 22.1 23.5 24.2 25.9 26.5 26.5 26.5 28.6]';
%! species = {'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' ...
%!            'ex' 'ex' 'ex' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' ...
%!            'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv'};
%!
%! stats = bootlm (pulse, {temp, species}, 'model', 'linear', ...
%!                           'continuous', 1, 'display', 'off', ...
%!                           'varnames', {'temp', 'species'}, 'seed', 1);
%!
%! assert (stats.CI_lower(2), 3.408042874444448, 1e-09);
%! assert (stats.CI_upper(2), 3.797462875271906, 1e-09);
%! assert (stats.CI_lower(3), -11.39708913283446, 1e-09);
%! assert (stats.CI_upper(3), -8.733493336263452, 1e-09);

%!test
%!
%! # Factorial design with continuous covariate. The data is from a study of the
%! # effects of treatment and exercise on stress reduction score after adjusting
%! # for age. Data from R datarium package).
%!
%! score = [95.6 82.2 97.2 96.4 81.4 83.6 89.4 83.8 83.3 85.7 ...
%!          97.2 78.2 78.9 91.8 86.9 84.1 88.6 89.8 87.3 85.4 ...
%!          81.8 65.8 68.1 70.0 69.9 75.1 72.3 70.9 71.5 72.5 ...
%!          84.9 96.1 94.6 82.5 90.7 87.0 86.8 93.3 87.6 92.4 ...
%!          100. 80.5 92.9 84.0 88.4 91.1 85.7 91.3 92.3 87.9 ...
%!          91.7 88.6 75.8 75.7 75.3 82.4 80.1 86.0 81.8 82.5]';
%! treatment = {'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' ...
%!              'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' ...
%!              'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' ...
%!              'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  ...
%!              'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  ...
%!              'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'}';
%! exercise = {'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  ...
%!             'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' ...
%!             'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  ...
%!             'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  ...
%!             'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' ...
%!             'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'}';
%! age = [59 65 70 66 61 65 57 61 58 55 62 61 60 59 55 57 60 63 62 57 ...
%!        58 56 57 59 59 60 55 53 55 58 68 62 61 54 59 63 60 67 60 67 ...
%!        75 54 57 62 65 60 58 61 65 57 56 58 58 58 52 53 60 62 61 61]';
%!
%! % ANOVA/ANCOVA statistics
%! [stats, bootstat, aovstat] = bootlm (score, {age, exercise, treatment}, ...
%!                            'model', [1 0 0; 0 1 0; 0 0 1; 0 1 1], ...
%!                            'continuous', 1, 'display', 'off', 'seed', 1, ...
%!                            'varnames', {'age', 'exercise', 'treatment'});
%!
%! assert (aovstat.PVAL(1), 0.0001, 1e-09);
%! assert (aovstat.PVAL(2), 0.0001, 1e-09);
%! assert (aovstat.PVAL(3), 0.00209853874900942, 1e-09);
%! assert (aovstat.PVAL(4), 0.0145576845409309, 1e-09);
%!
%! stats = bootlm (score, {age, exercise, treatment}, 'seed', 1, ...
%!                            'model', [1 0 0; 0 1 0; 0 0 1; 0 1 1], ...
%!                            'continuous', 1, 'display', 'off', ...
%!                            'varnames', {'age', 'exercise', 'treatment'});
%!
%! assert (stats.pval(6), 0.9605479032987221, 1e-09);
%! assert (stats.pval(7), 0.01418066878652798, 1e-09);
%! assert (stats.fpr(6), 0.5, 1e-09);
%! assert (stats.fpr(7), 0.1409314554632885, 1e-09);
%!
%! stats = bootlm (score, {age, exercise, treatment}, 'seed', 1, ...
%!                            'model', [1 0 0; 0 1 0; 0 0 1; 0 1 1], ...
%!                            'continuous', 1, 'display', 'off', ...
%!                            'varnames', {'age', 'exercise', 'treatment'}, ...
%!                            'dim', [2, 3]);
%!
%! assert (stats.estimate(1), 86.9787857062843,1e-09)
%! assert (stats.estimate(2), 86.9962428587431,1e-09)
%! assert (stats.estimate(3), 73.2754755236922,1e-09)
%! assert (stats.estimate(4), 88.5073652962921 ,1e-09)
%! assert (stats.estimate(5), 88.6798510137784,1e-09)
%! assert (stats.estimate(6), 83.02227960120982,1e-09)
%!
%! stats = bootlm (score, {age, exercise, treatment}, 'seed', 1, ...
%!                            'model', [1 0 0; 0 1 0; 0 0 1; 0 1 1], ...
%!                            'continuous', 1, 'display', 'off', ...
%!                            'varnames', {'age', 'exercise', 'treatment'}, ...
%!                            'dim', [2, 3], 'posthoc', 'trt_vs_ctrl');
%!
%! assert (stats.estimate(1), -0.0174571524588316,1e-09)
%! assert (stats.estimate(2), 13.7033101825921,1e-09)
%! assert (stats.estimate(3), -1.52857959000781,1e-09)
%! assert (stats.estimate(4), -1.70106530749405,1e-09)
%! assert (stats.estimate(5), 3.9565061050745,1e-09)

%!test
%!
%! # Unbalanced one-way design with custom, orthogonal contrasts.
%!
%! dv =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! g = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 ...
%!      4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%! C = [ 0.4001601  0.3333333  0.5  0.0
%!       0.4001601  0.3333333 -0.5  0.0
%!       0.4001601 -0.6666667  0.0  0.0
%!      -0.6002401  0.0000000  0.0  0.5
%!      -0.6002401  0.0000000  0.0 -0.5];
%!
%! stats = bootlm (dv, g, 'contrasts', C, 'varnames', 'score', 'seed', 1, ...
%!                          'alpha', 0.05, 'display', false);
%!
%! assert (stats.pval(2), 0.0001, 1e-09);
%! assert (stats.pval(3), 0.00189427105584975, 1e-09);
%! assert (stats.pval(4), 0.0001532987133628411, 1e-09);
%! assert (stats.pval(5), 0.0001, 1e-09);
%! assert (stats.fpr(2), 0.00249737757706675, 1e-09);
%! assert (stats.fpr(3), 0.03127029873554629, 1e-09);
%! assert (stats.fpr(4), 0.003646660191047087, 1e-09);
%! assert (stats.fpr(5), 0.00249737757706675, 1e-09);
%!
%! stats = bootlm (dv, g, 'contrasts', C, 'varnames', 'score', 'seed', 1, ...
%!                          'alpha', 0.05, 'display', false, 'dim', 1);
%!
%! assert (stats.CI_lower(1), 7.779565592818237, 1e-09);
%! assert (stats.CI_lower(2), 14.42536726599337, 1e-09);
%! assert (stats.CI_lower(3), 16.41718457146695, 1e-09);
%! assert (stats.CI_lower(4), 18.52263878670194, 1e-09);
%! assert (stats.CI_lower(5), 26.66171082767947, 1e-09);
%! assert (stats.CI_upper(1), 12.22043440718181, 1e-09);
%! assert (stats.CI_upper(2), 21.57463273400666, 1e-09);
%! assert (stats.CI_upper(3), 21.58281542853307, 1e-09);
%! assert (stats.CI_upper(4), 23.47764692758378, 1e-09);
%! assert (stats.CI_upper(5), 31.33851139454277, 1e-09);

%!test
%!
%! # One-way design.
%!
%! g = [1, 1, 1, 1, 1, 1, 1, 1, ...
%!      2, 2, 2, 2, 2, 2, 2, 2, ...
%!      3, 3, 3, 3, 3, 3, 3, 3]';
%! y = [13, 16, 16,  7, 11,  5,  1,  9, ...
%!      10, 25, 66, 43, 47, 56,  6, 39, ...
%!      11, 39, 26, 35, 25, 14, 24, 17]';
%!
%! stats = bootlm (y, g, 'display', false, 'dim', 1, 'posthoc', 'pairwise', ...
%!                       'seed', 1);
%!
%! assert (stats.pval(1), 0.02381212481394462, 1e-09);
%! assert (stats.pval(2), 0.009547350172112052, 1e-09);
%! assert (stats.pval(3), 0.1541408530918242, 1e-09);
%! assert (stats.fpr(1), 0.1254120362272493, 1e-09);
%! assert (stats.fpr(2), 0.04738586302975815, 1e-09);
%! assert (stats.fpr(3), 0.4392984660114189, 1e-09);
