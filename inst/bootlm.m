% -- statistics: bootlm (Y, GROUP)
% -- statistics: bootlm (Y, GROUP, NAME, VALUE)
% -- statistics: STATS = bootlm (...)
% -- statistics: [STATS, X] = bootlm (...)
%
%        Fits a linear model with categorical and/or continuous predictors (i.e.
%     independent variables) on a continuous outcome (i.e. dependent variable)
%     and computes the following statistics for each regression coefficient:
%          - name: the name(s) of the regression coefficient(s)
%          - coeff: the value of the regression coefficient(s)
%          - CI_lower: lower bound(s) of the 95% confidence interval (CI)
%          - CI_upper: upper bound(s) of the 95% confidence interval (CI)
%          - p-val: two-tailed p-value(s) for the parameter(s) being equal to 0
%        Confidence intervals and Null Hypothesis Significance Tests (NHSTs) for
%     the regression coefficients (H0 = 0) are calculated by wild bootstrap-t
%     and so are robust when normality and homoscedasticity cannot be assumed.
%     Note that the p-values are NOT adjusted for multiple comparisons.
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
%        • DIM is a scalar or vector specifying the dimension(s) over which
%          'bootlm' calculates and returns estimated marginal means instead of
%          regression coefficients. For example, the value [1 3] computes the
%          estimated marginal mean for each combination of the levels of the
%          first and third predictors. The default value is empty, which makes
%          'bootlm' return the statistics for for the model coefficients. If DIM
%          is, or includes, a continuous predictor then 'bootlm' will return an
%          error. The following statistics are returned when specifying 'dim':
%          - name: the name(s) of the estimated marginal mean(s)
%          - mean: the estimated marginal mean(s)
%          - CI_lower: lower bound(s) of the 95% confidence interval (CI)
%          - CI_upper: upper bound(s) of the 95% confidence interval (CI)
%          - n: the sample size used to estimate the mean
%
%     '[...] = bootlm (Y, GROUP, ..., 'continuous', CONTINUOUS)'
%
%        • CONTINUOUS is a vector of indices indicating which of the
%          columns (i.e. predictors) in GROUP should be treated as
%          continuous predictors rather than as categorical predictors.
%          The relationship between continuous predictors and the outcome
%          should be linear.
%
%     '[...] = bootlm (Y, GROUP, ..., 'model', MODELTYPE)'
%
%        • MODELTYPE can specified as one of the following:
%
%             • 'linear' (default): compute N main effects with no
%               interactions.
%
%             • 'interaction': compute N effects and N*(N-1) interactions
%
%             • 'full': compute the N main effects and interactions at
%               all levels
%
%             • a scalar integer: representing the maximum interaction
%               order
%
%             • a matrix of term definitions: each row is a term and
%               each column is a predictor
%
%               -- Example:
%               A two-way design with interaction would be: [1 0; 0 1; 1 1]
%
%     '[...] = bootlm (Y, GROUP, ..., 'varnames', VARNAMES)'
%
%        • VARNAMES must be a cell array of strings with each element
%          containing a predictor name for each column of GROUP. By default
%          (if not parsed as optional argument), VARNAMES are
%          'X1','X2','X3', etc.
%
%     '[...] = bootlm (Y, GROUP, ..., 'alpha', ALPHA)'
%
%        • ALPHA is numeric and sets the lower and upper bounds of the
%          confidence interval(s). The value(s) of ALPHA must be between
%          0 and 1. ALPHA can either be:
%
%             • scalar: To set the (nominal) central coverage of SYMMETRIC
%                  bootstrap-t confidence interval(s) to 100*(1-ALPHA)%.
%                  For example, 0.05 for a 95% confidence interval.
%
%             • vector: A pair of probabilities defining the (nominal) lower
%                  and upper bounds of ASYMMETRIC bootstrap-t confidence
%                  interval(s) as 100*(ALPHA(1))% and 100*(ALPHA(2))%
%                  respectively. For example, [.025, .975] for a 95%
%                  confidence interval.
%
%               The default value of ALPHA is the scalar: 0.05, for a symmetric
%               95% bootstrap-t confidence interval.
%
%     '[...] = bootlm (Y, GROUP, ..., 'display', DISPOPT)'
%
%        • DISPOPT can be either 'on' (or true, default) or 'off' (or false)
%          and controls the display of the model formula, a table of model
%          parameter estimates and a figure of diagnostic plots. The p-values
%          are formatted in APA-style.
%
%     '[...] = bootlm (Y, GROUP, ..., 'contrasts', CONTRASTS)'
%
%        • CONTRASTS can be specified as one of the following:
%
%             • A string corresponding to one of the built-in contrasts
%               listed below:
%
%                  • 'simple' or 'anova' (default): Simple (ANOVA) contrast
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
%                  • 'poly': Polynomial contrast coding for trend analysis.
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
%                  • 'helmert': Helmert contrast coding. The intercept
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
%                  • 'effect': Deviation effect coding. The intercept represents
%                    the grand mean. Each slope coefficient compares one level
%                    of a predictor (or interaction between predictors) with the
%                    grand mean. Note that a slope coefficient is omitted for
%                    the first level of the predictor(s) listed in the GROUP
%                    argument. The columns of this contrast coding scheme sum to
%                    zero. This type of contrast is ideal for nominal predictor
%                    variables when there is no obvious reference group.
%
%                  • 'sdif' or 'sdiff': Successive differences contrast coding.
%                    The intercept represents the grand mean. Each slope
%                    coefficient represents the difference between one level of
%                    a predictor (or interaction between predictors) to the
%                    previous one, where the order of the predictor levels is
%                    as they appear in the GROUP argument. The columns of this
%                    contrast coding coding scheme sum to zero. This type of
%                    contrast is ideal for ordinal predictor variables.
%
%                  • 'treatment': Treatment contrast (or dummy) coding. The
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
%             • A matrix containing a custom contrast coding scheme (i.e.
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
%        • Specifies the number of bootstrap resamples, where NBOOT must be a
%          positive integer. If empty, the default value of NBOOT is 10000.
%
%     '[...] = bootlm (Y, GROUP, ..., 'clustid', CLUSTID)'
%
%        • Specifies a vector or cell array of numbers or strings respectively
%          to be used as cluster labels or identifiers. Rows of the data with
%          the same CLUSTID value are treated as clusters with dependent errors.
%          If empty (default), no clustered resampling is performed and all
%          errors are treated as independent. The standard errors computed are
%          cluster robust.
%
%     '[...] = bootlm (Y, GROUP, ..., 'blocksz', BLOCKSZ)'
%
%        • Specifies a scalar, which sets the block size for bootstrapping when
%          the errors have serial dependence. Rows of the data within the same
%          block are treated as having dependent errors. If empty (default),
%          no clustered resampling is performed and all errors are treated
%          as independent. The standard errors computed are cluster robust.
%
%     '[...] = bootlm (Y, GROUP, ..., 'seed', SEED)' initialises the Mersenne
%     Twister random number generator using an integer SEED value so that
%     'bootlm' results are reproducible.
%
%     'bootlm' can return up to one output arguments:
%
%     'STATS = bootlm (...)' returns a structure with the following fields:
%     original, std_err, CI_lower, CI_upper, tstat, pval, fpr & name.
%
%     'STATS = bootlm (..., 'dim', DIM)' returns a structure with the following
%     fields: original, std_err, CI_lower, CI_upper, tstat, pval, fpr, name & n.
%
%     '[STATS, X] = bootlm (...)' also returns the design matrix for the linear 
%     model.
%
%     '[STATS, X, L] = bootlm (...)' also returns the hypothesis matrix used to
%     compute the estimated marginal means.
%
%  bootlm (version 2023.06.19)
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

function [STATS, X, L] = bootlm (Y, GROUP, varargin)

    if (nargin < 2)
      error (strcat (['bootlm usage: ''bootlm (Y, GROUP)''; '], ...
                      [' atleast 2 input arguments required']));
    end
    if (nargout > 3)
      error ('bootlm: Too many output arguments')
    end

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
    NBOOT = 10000;
    SEED = [];
    DEP = [];
    L = [];
    for idx = 3:2:nargin
      name = varargin{idx-2};
      value = varargin{idx-1};
      switch (lower (name))
        case 'model'
          MODELTYPE = value;
        case 'continuous'
          CONTINUOUS = value;
        case 'random'
          % RANDOM input argument is ignored
        case 'nested'
          error (strcat (['bootlm: NESTED not supported. Please use ''CLUSTID'''], ...
                         [' or ''BLOCKSZ'' input arguments instead.']));
        case 'sstype'
          % SSTYPE input argument is ignored
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
        case 'nboot'
          NBOOT = value;
        case 'seed'
          SEED = value;
        otherwise
          error (sprintf ('bootlm: parameter %s is not supported', name));
      end
    end

    % Evaluate continuous input argument
    if (isnumeric (CONTINUOUS))
      if (any (CONTINUOUS ~= abs (fix (CONTINUOUS))))
        error (strcat (['bootlm: the value provided for the CONTINUOUS'], ...
                       [' parameter must be a positive integer']));
      end
    else
      error (strcat (['bootlm: the value provided for the CONTINUOUS'], ...
                     [' parameter must be numeric']));
    end

    % Accomodate for different formats for GROUP
    % GROUP can be a matrix of numeric identifiers of a cell arrays
    % of strings or numeric identifiers
    N = size (GROUP, 2); % number of predictors
    n = numel (Y);       % total number of observations
    if (prod (size (Y)) ~= n)
      error ('bootlm: for ''bootlm (Y, GROUP)'', Y must be a vector');
    end
    if (numel (unique (CONTINUOUS)) > N)
      error (strcat (['bootlm: the number of predictors assigned as continuous'], ...
                 [' cannot exceed the number of predictors in GROUP']));
    end
    if (any ((CONTINUOUS > N) | any (CONTINUOUS <= 0)))
      error (strcat (['bootlm: one or more indices provided in the value'], ...
                 [' for the continuous parameter are out of range']));
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
              error ('bootlm: continuous predictors must be a numeric datatype');
            end
            tmp(:,j) = GROUP{j};
          end
        end
        GROUP = tmp;
      end
    end
    if (~ isempty (GROUP))
      if (size (GROUP,1) ~= n)
        error ('bootlm: GROUP must be a matrix with the same number of rows as Y');
      end
    end
    if (~ isempty (VARNAMES))
      if (iscell (VARNAMES))
        if (all (cellfun (@ischar, VARNAMES)))
          nvarnames = numel(VARNAMES);
        else
          error (strcat (['bootlm: all variable names must be character'], ...
                         [' or character arrays']));
        end
      elseif (ischar (VARNAMES))
        nvarnames = 1;
        VARNAMES = {VARNAMES};
      elseif (isstring (VARNAMES))
        nvarnames = 1;
        VARNAMES = {char(VARNAMES)};
      else
        error (strcat (['bootlm: varnames is not of a valid type. Must be a cell'], ...
               [' array of character arrays, character array or string']));
      end
    else
      nvarnames = N;
      VARNAMES = arrayfun(@(x) ['X',num2str(x)], 1:N, 'UniformOutput', 0);
    end
    if (nvarnames ~= N)
      error (strcat (['bootlm: number of variable names is not equal'], ...
                     [' to the number of grouping variables']));
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
              'Note that the CONTRASTS for predictor %u do not sum to zero', i));
            end
            % Check whether contrasts are orthogonal
            if (any (abs (reshape (corr (CONTRASTS{i}) - ...
                                          eye (size (CONTRASTS{i}, 2)), [], 1))...
                                        > eps ('single')))
              warning (sprintf ( ...
              'Note that the CONTRASTS for predictor %u are not orthogonal', i));
            end
          else
            if (~ ismember (lower (CONTRASTS{i}), ...
                            {'simple','anova','poly','helmert','effect',...
                              'sdif','sdiff','treatment'}))
              error (strcat(['bootlm: valid built-in contrasts are:'], ...
                            [' ''simple'', ''poly'', ''helmert'','],...
                            ['''effect'', ''sdif'' or ''treatment''']));
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
    msg = strcat (['bootlm: the number of columns in the term definitions'], ...
                  [' cannot exceed the number of columns of GROUP']);
    if (ischar (MODELTYPE))
      switch (lower (MODELTYPE))
        case 'linear'
          MODELTYPE = 1;
        case {'interaction','interactions'}
          MODELTYPE = 2;
        case 'full'
          MODELTYPE = N;
        otherwise
          error ('bootlm: model type not recognised');
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
        error (strcat (['bootlm: elements of the model terms matrix'], ...
                       [' must be either 0 or 1']));
      end
      TERMS = logical (MODELTYPE);
    end
    % Evaluate terms matrix
    Ng = sum (TERMS, 2);
    if (any (diff (Ng) < 0))
      error (strcat (['bootlm: the model terms matrix must list main'], ...
                     [' effects above/before interactions']));
    end
    % Evaluate terms
    Nm = sum (Ng == 1);
    Nx = sum (Ng > 1);
    Nt = Nm + Nx;
    if (any (any (TERMS(1:Nm,:), 1) ~= any (TERMS, 1)))
      error (strcat (['bootlm: all predictors involved in interactions'], ...
                     [' must have a main effect']));
    end

    % Create design matrix
    [X, grpnames, nlevels, df, termcols, coeffnames, vmeans, gid, CONTRASTS, center_continuous] = ...
      mDesignMatrix (GROUP, TERMS, CONTINUOUS, CONTRASTS, VARNAMES, n, Nm, Nx, Ng, cont_vec);
    dft = n - 1;
    dfe = dft - sum (df);

    % If applicable, create hypothesis matrix, names and compute sample sizes
    if (~ isempty (DIM))
      H = X;
      for i = 1:Nt
        if ( all (DIM ~= i) )
          H{i+1}(:,:) = 0;
        end
      end
      H = cell2mat (H);
      L = unique (H, 'rows', 'stable').';
    end

    % Fit linear model, and calculate sums-of-squares
    X = cell2mat (X);
    [b, sse, resid, hat, ucov] = lmfit (X, Y);

    % Prepare model formula 
    formula = sprintf ('Y ~ 1');  % Initialise model formula
    for i = 1:Nt
      str = sprintf ('%s*', VARNAMES{find (TERMS(i,:))});
      % Append model term to formula
      str = regexprep (str, '\\*', ':');
      % Fixed effect term
      formula = sprintf ('%s + %s', formula, str(1:end-1));
    end

    % Use bootstrap methods to calculate statistics
    if isempty (DIM)

      % Model coefficients
      STATS = bootwild (Y, X, DEP, NBOOT, ALPHA, SEED);

      % Assign names for model coefficients
      NAMES = vertcat (coeffnames{:});
      STATS.name = NAMES;

    else

      % Model estimated marginal means
      STATS = bootwild (Y, X, DEP, NBOOT, ALPHA, SEED, L);

      % Create names for estimated marginal means
      idx = cellfun (@(l) find (all (bsxfun (@eq, H, l), 2), 1), num2cell (L', 2));
      Ne = size (L, 2);
      Np = numel (DIM);
      NAMES = cell (Ne, 1);
      for i = 1 : Ne
        str = '';
        for j = 1 : Np
          str = sprintf('%s%s=%s, ', str, ...
                    num2str (VARNAMES{DIM(j)}), ...
                    num2str (grpnames{DIM(j)}{gid(idx(i),DIM(j))}));
        end
        NAMES{i} = str(1:end-2);
        str = '';
      end
      STATS.name = NAMES;

      % Compute sample sizes and add them to the output structure
      U = unique (gid(:,DIM), 'rows', 'stable');
      STATS.n = cellfun (@(u) sum (all (gid(:,DIM) == u, 2)), num2cell (U, 2));

    end

    % Print table of model coefficients and make figure of diagnostic plots
    switch (lower (DISPLAY))

      case {'on', true}

        % Print model formula
        fprintf('\nMODEL FORMULA (based on Wilkinson''s notation):\n\n%s\n', formula);

        % If applicable, print parameter estimates (a.k.a contrasts) for fixed effects
        % Parameter estimates correspond to the contrasts we set.
        if (isempty (DIM))
          fprintf('\nMODEL COEFFICIENTS\n\n');
          fprintf('name                                   coeff       CI_lower    CI_upper    p-val\n');
          fprintf('--------------------------------------------------------------------------------\n');
        else
          fprintf('\nMODEL ESTIMATED MARGINAL MEANS\n\n');
          fprintf('name                                   mean        CI_lower    CI_upper        n\n');
          fprintf('--------------------------------------------------------------------------------\n');
        end
        for j = 1:size (NAMES, 1)
          if (isempty (DIM))
            fprintf ('%-37s  %#-+10.4g  %#-+10.4g  %#-+10.4g', ...
                     NAMES{j}(1:min(end,37)), STATS.original(j), STATS.CI_lower(j), STATS.CI_upper(j));
            if (STATS.pval(j) <= 0.001)
              fprintf ('  <.001\n');
            elseif (STATS.pval(j) < 0.9995)
              fprintf ('   .%03u\n', round (STATS.pval(j) * 1e+03));
            elseif (isnan (STATS.pval(j)))
              fprintf ('    NaN\n');
            else
              fprintf ('   1.000\n');
            end
          else
            fprintf ('%-37s  %#-+10.4g  %#-+10.4g  %#-+10.4g  %5u\n', ...
                     NAMES{j}(1:min(end,37)), STATS.original(j), STATS.CI_lower(j), STATS.CI_upper(j), STATS.n(j));
          end
        end
        fprintf('\n');

        % Make figure of diagnostic plots
        figure("Name", "Diagnostic plots: Standardized Model Residuals");
        resid = resid;
        std_resid = (resid - mean (resid)) ./ std (resid, 1);
        fit = X * b;
      
        % Normal probability plot
        subplot (2, 2, 1);
        normplot (std_resid);
        xlabel ("Standardized Residuals");
      
        % Checks for homoskedasticity assumption
        subplot (2, 2, 2);
        plot (fit, std_resid, "b+");
        xlabel ("Fitted values");
        ylabel ("Standardized Residuals");
        title ("Residuals vs Fitted Values")
        ax1_xlim = get (gca, 'XLim');
        hold on; plot (ax1_xlim, zeros (1, 2), "r--"); grid ("on"); hold off;
      
        % Checks for outliers and heteroskedasticity
        subplot (2, 2, 3);
        plot (fit, sqrt (abs (std_resid)), "b+");
        xlabel ("Fitted values");
        ylabel ("sqrt ( |Standardized Residuals| )");
        title ("Spread-Location Plot")
        ax2_xlim = get (gca, 'XLim');
        hold on; 
        plot (ax2_xlim, ones (1, 2) * sqrt (2), "b--");
        plot (ax2_xlim, ones (1, 2) * sqrt (3), "m--"); 
        plot (ax2_xlim, ones (1, 2) * sqrt (4), "r--");
        hold off;
      
        % Check for influential outliers
        subplot (2, 2, 4);
        plot (diag (hat), std_resid, "b+");
        xlabel ("Leverage")
        ylabel ("Standardized Residuals");
        title ("Influential values")
        ax3_xlim = get (gca, 'XLim');
        ax3_ylim = get (gca, 'YLim');
        hold on; plot (ax3_xlim, zeros (1, 2), "r--"); grid ("on"); hold off;
        ylim (ax3_ylim);
        set (findall ( gcf, '-property', 'FontSize'), 'FontSize', 7)

      case {'off', false}

        % do nothing

      otherwise

        error ('bootlm: wrong value for ''display'' parameter.');

  end

end

function [X, levels, nlevels, df, termcols, coeffnames, vmeans, gid, CONTRASTS, center_continuous] = ...
         mDesignMatrix (GROUP, TERMS, CONTINUOUS, CONTRASTS, VARNAMES, n, Nm, Nx, Ng, cont_vec)

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
      levels{j} = unique (GROUP(:,j), 'stable');
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
          error (strcat (['bootlm: the number of rows in the contrast'], ...
                       [' matrices should equal the number of predictor levels']));
        end
        if (~ all (size (CONTRASTS{j},2) == df(j)))
          error (strcat (['bootlm: the number of columns in each contrast'], ...
                  [' matrix should equal the degrees of freedom (i.e.'], ...
                  [' number of levels minus 1) for that predictor']));
        end
        if (~ all (any (CONTRASTS{j})))
          error (strcat (['bootlm: a contrast must be coded in each'], ...
                         [' column of the contrast matrices']));
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
      X{1+j} = cell2mat (cellfun (func, num2cell (C, 1), 'UniformOutput', false));
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

% BUILT IN CONTRAST CODING FUNCTIONS

function C = contr_simple (N)

  % Create contrast matrix (of doubles) using simple (ANOVA) contrast coding
  % These contrasts are centered (i.e. sum to 0)
  % Ideal for unordered predictors, with comparison to a reference level
  % The first predictor level is the reference level
  C =  cat (1, zeros (1,N-1), eye(N-1)) - 1/N;

end

function C = contr_poly (N)

  % Create contrast matrix (of doubles) using polynomial contrast coding
  % for trend analysis of ordered categorical predictor levels
  % These contrasts are orthogonal and centered (i.e. sum to 0)
  % Ideal for ordered predictors
  [C, jnk] = qr (bsxfun (@power, [1:N]' - mean ([1:N]'), [0:N-1]));
  C(:,1) = [];
  s = ones (1, N-1);
  s(1:2:N-1) = s(1:2:N-1) * -1;
  f = (sign(C(1,:)) ~= s);
  C(:,f) = C(:,f) * -1;

end

function C = contr_helmert (N)

  % Create contrast matrix (of doubles) using Helmert coding contrasts
  % These contrasts are orthogonal and centered (i.e. sum to 0)
  C = cat (1, tril (-ones (N-1), -1) + diag (N-1:-1:1), ...
              -ones (1, N-1)) ./ (N:-1:2);

end

function C = contr_sum (N)

  % Create contrast matrix (of doubles) using deviation effect coding
  % These contrasts are centered (i.e. sum to 0)
  C =  cat (1, - (ones (1,N-1)), eye (N-1));

end

function C = contr_sdif (N)

  % Create contrast matrix (of doubles) using successive differences coding
  % These contrasts are centered (i.e. sum to 0)
  C =  tril (ones (N, N - 1), -1) - ones (N, 1) / N * [N - 1 : -1 : 1];

end

function C = contr_treatment (N)

  % Create contrast matrix (of doubles) using treatment contrast coding
  % Ideal for unordered predictors, with comparison to a reference level
  % The first predictor level is the reference level
  C =  cat (1, zeros (1,N-1), eye(N-1));

end


% FUNCTION TO FIT THE LINEAR MODEL

function [b, sse, resid, hat, ucov] = lmfit (X, Y)

  % Get model coefficients by solving the linear equation by QR decomposition
  % The number of free parameters (i.e. intercept + coefficients) is equal
  % to n - dfe (i.e. the number of columns in X).
  [Q, R] = qr (X, 0);
  b = R \ Q' * Y;            % Instead of pinv (X' * X) * (X' * y);

  % Calculate the Hat matrix (i.e. X*(X'*X)^−1*X')
  hat = Q * Q';

  % Get fitted values
  fit = hat * Y;             % Instead of X * b;

  % Get residuals from the fit
  resid = Y - fit;

  % Calculate the residual sums-of-squares
  sse = sum (resid.^2);

  % Calculate the unscaled covariance matrix (i.e. inv (X'*X ))
  if (nargout > 4)
    ucov = R \ Q' / X';
  end

end


%!demo
%!
%! # Two-sample unpaired test on independent samples (equivalent to Student's
%! # t-test). Note that the absolute value of t-statistic can be obtained by
%! # taking the square root of the reported F statistic. In this example,
%! # t = sqrt (1.44) = 1.20.
%!
%! score = [54 23 45 54 45 43 34 65 77 46 65]';
%! gender = {'male' 'male' 'male' 'male' 'male' 'female' 'female' 'female' ...
%!           'female' 'female' 'female'}';
%!
%! STATS = bootlm (score, gender, 'display', 'on', 'varnames', 'gender');

%!demo
%!
%! # Two-sample paired test on dependent or matched samples equivalent to a
%! # paired t-test. As for the first example, the t-statistic can be obtained by
%! # taking the square root of the reported F statistic.
%!
%! score = [4.5 5.6; 3.7 6.4; 5.3 6.4; 5.4 6.0; 3.9 5.7]';
%! treatment = {'before' 'after'; 'before' 'after'; 'before' 'after';
%!              'before' 'after'; 'before' 'after'}';
%! subject = {'GS' 'GS'; 'JM' 'JM'; 'HM' 'HM'; 'JW' 'JW'; 'PS' 'PS'}';
%!
%! STATS = bootlm (score(:), {treatment(:), subject(:)}, ...
%!                            'model', 'linear', 'display', 'on', ...
%!                            'varnames', {'treatment', 'subject'});

%!demo
%!
%! # One-way design on the data from a study on the strength of structural beams,
%! # in Hogg and Ledolter (1987) Engineering Statistics. New York: MacMillan
%!
%! strength = [82 86 79 83 84 85 86 87 74 82 ...
%!            78 75 76 77 79 79 77 78 82 79]';
%! alloy = {'st','st','st','st','st','st','st','st', ...
%!          'al1','al1','al1','al1','al1','al1', ...
%!          'al2','al2','al2','al2','al2','al2'}';
%!
%! STATS = bootlm (strength, alloy, 'display', 'on', 'varnames', 'alloy');

%!demo
%!
%! # One-way repeated measures design on the data from a study on the number of
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
%! STATS = bootlm (words(:), {seconds(:), subject(:)}, ...
%!                            'model', 'linear', 'display', 'on', ...
%!                            'varnames', {'seconds', 'subject'});

%!demo
%!
%! # Balanced two-way design with interaction on the data from a study of popcorn
%! # brands and popper types, in Hogg and Ledolter (1987) Engineering Statistics.
%! # New York: MacMillan
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
%! STATS = bootlm (popcorn(:), {brands(:), popper(:)}, ...
%!                            'display', 'on', 'model', 'full', ...
%!                            'varnames', {'brands', 'popper'});

%!demo
%!
%! # Unbalanced two-way design (2x2) on the data from a study on the effects of
%! # gender and having a college degree on salaries of company employees,
%! # in Maxwell, Delaney and Kelly (2018): Chapter 7, Table 15
%!
%! salary = [24 26 25 24 27 24 27 23 15 17 20 16, ...
%!           25 29 27 19 18 21 20 21 22 19]';
%! gender = {'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f'...
%!           'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm'}';
%! degree = [1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 0 0 0 0 0 0 0]';
%!
%! STATS = bootlm (salary, {gender, degree}, 'model', 'full', ...
%!                            'display', 'on', 'varnames', ...
%!                            {'gender', 'degree'});

%!demo
%!
%! # Unbalanced two-way design (3x2) on the data from a study of the effect of
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
%! STATS = bootlm (babble, {sugar, milk}, 'model', 'full', 'display', 'on', ...
%!                                        'varnames', {'sugar', 'milk'});

%!demo
%!
%! # Unbalanced three-way design (3x2x2) on the data from a study of the effects
%! # of three different drugs, biofeedback and diet on patient blood pressure,
%! # adapted* from Maxwell, Delaney and Kelly (2018): Chapter 8, Table 12
%! # * Missing values introduced to make the sample sizes unequal to test the
%! #   calculation of different types of sums-of-squares
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
%!       189 194 217 206 199 195 171 173 196 199 180 NaN;
%!       180 187 199 170 204 194 162 184 183 156 180 173 ...
%!       202 228 190 206 224 204 205 199 170 160 NaN NaN];
%!
%! STATS = bootlm (BP(:), {drug(:), feedback(:), diet(:)}, ...
%!                                    'model', 'full', ...
%!                                    'display', 'on', ...
%!                                    'varnames', {'drug', 'feedback', 'diet'});

%!demo
%!
%! # Balanced three-way design (2x2x2) with one of the predictors being a
%! # blocking factor. The data is from a randomized block design study on the
%! # effects of antioxidant treatment on glutathione-S-transferase (GST) levels
%! # in different mouse strains, from Festing (2014), ILAR Journal 55(3):427-476.
%!
%! measurement = [444 614 423 625 408  856 447 719 ...
%!                764 831 586 782 609 1002 606 766]';
%! strain= {'NIH','NIH','BALB/C','BALB/C','A/J','A/J','129/Ola','129/Ola', ...
%!          'NIH','NIH','BALB/C','BALB/C','A/J','A/J','129/Ola','129/Ola'}';
%! treatment={'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T'}';
%! block = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2]';
%!
%! STATS = bootlm (measurement/10, {strain, treatment, block}, ...
%!                            'model', [1 0 0; 0 1 0; 0 0 1; 1 1 0], ...
%!                            'varnames', {'strain', 'treatment', 'block'}, ...
%!                            'display', 'on');

%!demo
%!
%! # One-way design with continuous covariate on data from a study of the
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
%! STATS = bootlm (pulse, {species, temp}, 'model', 'linear', ...
%!                           'continuous', 2, 'display', 'on', ...
%!                           'varnames', {'species', 'temp'});

%!demo
%!
%! # Factorial design with continuous covariate on data from a study of the
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
%! STATS = bootlm (score, {treatment, exercise, age}, ...
%!                            'model', [1 0 0; 0 1 0; 0 0 1; 1 1 0], ...
%!                            'continuous', 3, 'display', 'on', ...
%!                            'varnames', {'treatment', 'exercise', 'age'});

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
%! STATS = bootlm (dv, g, 'contrasts', C, 'varnames', 'score', ...
%!                          'alpha', 0.05, 'display', true);
%!
%! STATS = bootlm (dv, g, 'contrasts', C, 'varnames', 'score', ...
%!                          'alpha', 0.05, 'display', true, 'dim', 1);

%!test
%!
%! # Two-sample unpaired test on independent samples (equivalent to Student's
%! # t-test). Note that the absolute value of t-statistic can be obtained by
%! # taking the square root of the reported F statistic. In this example,
%! # t = sqrt (1.44) = 1.20.
%!
%! score = [54 23 45 54 45 43 34 65 77 46 65]';
%! gender = {'male' 'male' 'male' 'male' 'male' 'female' 'female' 'female' ...
%!           'female' 'female' 'female'}';
%!
%! STATS = bootlm (score, gender, 'display', 'off', 'varnames', 'gender');

%!test
%!
%! # Two-sample paired test on dependent or matched samples equivalent to a
%! # paired t-test. As for the first example, the t-statistic can be obtained by
%! # taking the square root of the reported F statistic.
%!
%! score = [4.5 5.6; 3.7 6.4; 5.3 6.4; 5.4 6.0; 3.9 5.7]';
%! treatment = {'before' 'after'; 'before' 'after'; 'before' 'after';
%!              'before' 'after'; 'before' 'after'}';
%! subject = {'GS' 'GS'; 'JM' 'JM'; 'HM' 'HM'; 'JW' 'JW'; 'PS' 'PS'}';
%!
%! STATS = bootlm (score(:), {treatment(:), subject(:)}, ...
%!                            'model', 'linear', 'display', 'off', ...
%!                            'varnames', {'treatment', 'subject'});

%!test
%!
%! # One-way design on the data from a study on the strength of structural beams,
%! # in Hogg and Ledolter (1987) Engineering Statistics. New York: MacMillan
%!
%! strength = [82 86 79 83 84 85 86 87 74 82 ...
%!            78 75 76 77 79 79 77 78 82 79]';
%! alloy = {'st','st','st','st','st','st','st','st', ...
%!          'al1','al1','al1','al1','al1','al1', ...
%!          'al2','al2','al2','al2','al2','al2'}';
%!
%! STATS = bootlm (strength, alloy, 'display', 'off', 'varnames', 'alloy');

%!test
%!
%! # One-way repeated measures design on the data from a study on the number of
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
%! STATS = bootlm (words(:), {seconds(:), subject(:)}, ...
%!                            'model', 'linear', 'display', 'off', ...
%!                            'varnames', {'seconds', 'subject'});

%!test
%!
%! # Balanced two-way design with interaction on the data from a study of popcorn
%! # brands and popper types, in Hogg and Ledolter (1987) Engineering Statistics.
%! # New York: MacMillan
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
%! STATS = bootlm (popcorn(:), {brands(:), popper(:)}, ...
%!                            'display', 'off', 'model', 'full', ...
%!                            'varnames', {'brands', 'popper'});

%!test
%!
%! # Unbalanced two-way design (2x2) on the data from a study on the effects of
%! # gender and having a college degree on salaries of company employees,
%! # in Maxwell, Delaney and Kelly (2018): Chapter 7, Table 15
%!
%! salary = [24 26 25 24 27 24 27 23 15 17 20 16, ...
%!           25 29 27 19 18 21 20 21 22 19]';
%! gender = {'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f'...
%!           'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm'}';
%! degree = [1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 0 0 0 0 0 0 0]';
%!
%! STATS = bootlm (salary, {gender, degree}, 'model', 'full', ...
%!                            'display', 'off', 'varnames', ...
%!                            {'gender', 'degree'});

%!test
%!
%! # Unbalanced two-way design (3x2) on the data from a study of the effect of
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
%! STATS = bootlm (babble, {sugar, milk}, 'model', 'full', 'display', 'off', ...
%!                                        'varnames', {'sugar', 'milk'});

%!test
%!
%! # Unbalanced three-way design (3x2x2) on the data from a study of the effects
%! # of three different drugs, biofeedback and diet on patient blood pressure,
%! # adapted* from Maxwell, Delaney and Kelly (2018): Chapter 8, Table 12
%! # * Missing values introduced to make the sample sizes unequal to test the
%! #   calculation of different types of sums-of-squares
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
%!       189 194 217 206 199 195 171 173 196 199 180 NaN;
%!       180 187 199 170 204 194 162 184 183 156 180 173 ...
%!       202 228 190 206 224 204 205 199 170 160 NaN NaN];
%!
%! STATS = bootlm (BP(:), {drug(:), feedback(:), diet(:)}, ...
%!                                    'model', 'full', ...
%!                                    'display', 'off', ...
%!                                    'varnames', {'drug', 'feedback', 'diet'});

%!test
%!
%! # Balanced three-way design (2x2x2) with one of the predictors being a
%! # blocking factor. The data is from a randomized block design study on the
%! # effects of antioxidant treatment on glutathione-S-transferase (GST) levels
%! # in different mouse strains, from Festing (2014), ILAR Journal 55(3):427-476.
%!
%! measurement = [444 614 423 625 408  856 447 719 ...
%!                764 831 586 782 609 1002 606 766]';
%! strain= {'NIH','NIH','BALB/C','BALB/C','A/J','A/J','129/Ola','129/Ola', ...
%!          'NIH','NIH','BALB/C','BALB/C','A/J','A/J','129/Ola','129/Ola'}';
%! treatment={'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T'}';
%! block = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2]';
%!
%! STATS = bootlm (measurement/10, {strain, treatment, block}, ...
%!                            'model', [1 0 0; 0 1 0; 0 0 1; 1 1 0], ...
%!                            'varnames', {'strain', 'treatment', 'block'}, ...
%!                            'display', 'off');

%!test
%!
%! # One-way design with continuous covariate on data from a study of the
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
%! STATS = bootlm (pulse, {species, temp}, 'model', 'linear', ...
%!                           'continuous', 2, 'display', 'off', ...
%!                           'varnames', {'species', 'temp'});

%!test
%!
%! # Factorial design with continuous covariate on data from a study of the
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
%! STATS = bootlm (score, {treatment, exercise, age}, ...
%!                            'model', [1 0 0; 0 1 0; 0 0 1; 1 1 0], ...
%!                            'continuous', 3, 'display', 'off', ...
%!                            'varnames', {'treatment', 'exercise', 'age'});

%!test
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
%! STATS = bootlm (dv, g, 'contrasts', C, 'varnames', 'score', ...
%!                          'alpha', 0.05, 'display', false);
%!
%! STATS = bootlm (dv, g, 'contrasts', C, 'varnames', 'score', ...
%!                          'alpha', 0.05, 'display', false, 'dim', 1);
