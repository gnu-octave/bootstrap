% -- Function File: bootnhst (DATA, GROUP)
% -- Function File: bootnhst (..., NAME, VALUE)
% -- Function File: bootnhst (..., 'bootfun', BOOTFUN)
% -- Function File: bootnhst (..., 'nboot', NBOOT)
% -- Function File: bootnhst (..., 'ref', REF)
% -- Function File: bootnhst (..., 'alpha', ALPHA)
% -- Function File: bootnhst (..., 'Options', PAROPT)
% -- Function File: PVAL = bootnhst (DATA, GROUP, ...)
% -- Function File: [PVAL, C] = bootnhst (DATA, GROUP, ...)
% -- Function File: [PVAL, C, STATS] = bootnhst (DATA, GROUP, ...)
% -- Function File: [...] = bootnhst (..., 'display', DISPLAYOPT)
%
%     'bootnhst (DATA, GROUP)' performs a bootstrap version of a randomization
%     test [1] for comparing independent samples of data in a one-way layout,
%     where the categorization of the DATA values are defined by the labels in
%     GROUP. DATA must be a numeric column vector or matrix, and GROUP must be
%     a vector or cell array with the same number of rows as DATA. Pairwise
%     posthoc tests are automatically computed by the single-step maximum t-
%     statistic (maxT) procedure, which controls the family-wise error rate
%     (FWER) in a manner analagous to the Tukey-Kramer Honest Significance
%     Difference test. The omnibus test represents the smallest of the
%     multiplicity-adjusted p-values. The results are displayed as a pretty
%     table and the differences between groups are plotted along with the
%     symmetic 95% bootstrap-t confidence intervals (CI). The colours of the
%     markers and error bars depend on the value of the multiplicity-adjusted
%     p-values: red if p < .05, or blue if p > .05. All of the p-values reported
%     represent the outcome of two-tailed tests. 
%
%     bootnhst can take a number of optional parameters as NAME-VALUE pairs:
%
%     'bootnhst (..., 'bootfun', BOOTFUN)' also specifies BOOTFUN: the function
%     calculated on the original sample and the bootstrap resamples. BOOTFUN
%     must be either a:
%        o function handle or anonymous function,
%        o string of function name, or
%        o a cell array where the first cell is one of the above function
%          definitions and the remaining cells are (additional) input arguments 
%          to that function (other than the data arguments).
%        In all cases, BOOTFUN must take DATA for the initial input argument(s).
%        BOOTFUN must calculate a statistic representative of the finite data
%        sample; it should NOT be an estimate of a population parameter (unless
%        they are one of the same). By default, BOOTFUN is @mean. If a robust
%        alternative to the mean is required, set BOOTFUN to 'robust' to
%        implement a smoothed version of the median (a.k.a. @smoothmedian). 
%
%     'bootnhst (..., 'nboot', NBOOT)' is a scalar or a vector of upto two
%     positive integers indicating the number of resamples for the first
%     (bootstrap) and second (bootknife) levels of iterated resampling. If NBOOT
%     is a scalar value, or if NBOOT(2) is set to 0, then standard errors are
%     calculated either without resampling (if BOOTFUN @mean) or using Tukey's
%     jackknife. This implementation of jackknife requires the Statistics
%     package/toolbox. The default value of NBOOT is the vector: [999,99].
%
%     'bootnhst (..., 'ref', REF)' sets the GROUP to use as the reference group
%     for the post hoc tests. If REF is a recognised member of GROUP, then the
%     maxT procedure for treatment versus reference controls the family-wise
%     error rate (FWER) in a manner analagous to Dunnett's post hoc tests.
%
%     'bootnhst (..., 'alpha', ALPHA)' specifies the two-tailed significance
%     level for CI coverage. The default value of ALPHA is 0.05 for 95%
%     confidence intervals.
%
%     'bootnhst (..., 'Options', PAROPT)' specifies options that govern if
%     and how to perform bootstrap iterations using multiple processors (if the
%     Parallel Computing Toolbox or Octave Parallel package is available). This
%     argument is a structure with the following recognised fields:
%        o 'UseParallel':  If true, use parallel processes to accelerate
%                          bootstrap computations on multicore machines,
%                          specifically non-vectorized function evaluations,
%                          double bootstrap resampling and jackknife function
%                          evaluations. Default is false for serial computation.
%                          In MATLAB, the default is true if a parallel pool
%                          has already been started. 
%        o 'nproc':        nproc sets the number of parallel processes
%
%     'PVAL = bootnhst (DATA, GROUP, ...)' returns the p-value for the omnibus
%     hypothesis test. Note that the p-value returned will be truncated at the
%     resolution limit determined by the number of bootstrap replicates,
%     specifically 1/NBOOT(1). 
%
%     '[PVAL, C] = bootnhst (DATA, GROUP, ...)' also returns a 9 column matrix
%     that summarises post hoc test results. The columns of C are:
%       - column 1: reference GROUP number
%       - column 2: test GROUP number
%       - column 3: value of BOOTFUN evaluated for the reference GROUP
%       - column 4: value of BOOTFUN evaluated for the test GROUP
%       - column 5: the difference between the groups (columns 4 minus column 3)
%       - column 6: t-ratio
%       - column 7: multiplicity-adjusted p-value
%       - column 8: LOWER bound of the 100*(1-ALPHA)% bootstrap-t CI
%       - column 9: UPPER bound of the 100*(1-ALPHA)% bootstrap-t CI
%
%     '[PVAL, C, STATS] = bootnhst (DATA, GROUP, ...)' also returns a structure 
%     containing additional statistics. The stats structure contains the 
%     following fields:
%
%       gnames   - group names used in the GROUP input argument. The index of 
%                  gnames corresponds to the numbers used to identify GROUPs
%                  in columns 1 and 2 of the output argument c
%       ref      - index of the reference group
%       groups   - group index and BOOTFUN for each group with sample size,
%                  standard error and CI that start to overlap at a multiplicity
%                  adjusted p-value of approximately 0.05
%       Var      - weighted mean (pooled) sampling variance
%       maxT     - omnibus test statistic (maxT) 
%       df       - degrees of freedom (if bootfun is 'mean' or 'robust')
%       nboot    - number of bootstrap resamples
%       alpha    - two-tailed significance level for the CI reported in c.
%       bootstat - test statistic computed for each bootstrap resample 
%
%     '[...] = bootnhst (..., 'display', DISPLAYOPT)' a logical value (true 
%      or false) used to specify whether to display the results and plot the
%      graph in addition to creating the output arguments. The default is true.
%
%  BIBLIOGRAPHY:
%  [1] Efron, and Tibshirani (1993) An Introduction to the Bootstrap. 
%        New York, NY: Chapman & Hall
%
%  bootnhst (version 2023.07.04)
%  Bootstrap Null Hypothesis Significance Test
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
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.


function [pval, c, stats] = bootnhst (data, group, varargin)

  % Evaluate the number of function arguments
  if (nargin < 2)
    error (cat (2, 'bootnhst usage: ''bootnhst (DATA, GROUP, varargin)'';', ...
                   ' atleast 2 input arguments required'))
  end

  % Store local functions in a stucture for parallel processes
  localfunc = struct ('maxstat',@maxstat, ...
                      'bootcdf',@bootcdf);

  % Check if running in Octave (else assume Matlab)
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));
  
  % Apply defaults
  bootfun = 'mean';
  nboot = [999,99];
  ref = [];
  alpha = 0.05;
  DisplayOpt = true;
  paropt = struct;
  paropt.UseParallel = false;
  if (ISOCTAVE)
    paropt.nproc = nproc;
  else
    paropt.nproc = feature ('numcores');
  end

  % Fetch extra input arguments
  argin3 = varargin;
  narg = numel (argin3);
  if (narg > 1)
    while ischar (argin3{end-1})
      if (strcmpi (argin3{end-1},'bootfun'))
        bootfun = argin3{end};
      elseif (strcmpi (argin3{end-1},'nboot'))
        nboot = argin3{end};
      elseif (strcmpi (argin3{end-1},'ref'))
        ref = argin3{end};
      elseif (any (strcmpi (argin3{end-1},{'Options','Option'})))
        paropt = argin3{end};
      elseif (strcmpi (argin3{end-1},'alpha'))
        alpha = argin3{end};
      elseif (any (strcmpi (argin3{end-1},{'DisplayOpt','Display'})))
        DisplayOpt = argin3{end};
      else
        error ('bootnhst: Unrecognised input argument to bootnhst')
      end
      argin3 = {argin3{1:end-2}};
      narg = numel (argin3);
      if (narg < 1)
        break
      end
    end
  end

  % Error checking
  % Check and process bootnhst input arguments
  nvar = size (data,2);
  if (nargin < 2)
    error ('bootnhst: Requires atleast two input arguments');
  end
  if (ischar (group))
    group = cellstr (group);
  end
  if ((size (group, 1)>1) && (size (data, 1) ~= size (group, 1)))
    error ('bootnhst: DATA and GROUP must have the same number of rows')
  end
  if (iscell (group))
    if (~ iscellstr (group))
      group = cell2mat (group);
    end
  end
  if (iscell (bootfun))
    if (ischar (bootfun{1}))
      % Convert character string of a function name to a function handle
      bootfun_str = bootfun{1};
      func = str2func (bootfun{1});
    else
      bootfun_str = func2str (bootfun{1});
      func = bootfun{1};
    end
    args = bootfun(2:end);
    bootfun = @(varargin) func (varargin{:}, args{:});
  elseif (ischar (bootfun))
    % Convert character string of a function name to a function handle
    bootfun_str = bootfun;
    if (strcmpi (bootfun_str, 'robust'))
      bootfun = @smoothmedian;
    else
      bootfun = str2func (bootfun);
    end
  elseif (isa (bootfun, 'function_handle'))
    bootfun_str = func2str (bootfun);
  else
    error ('bootnhst: BOOTFUN must be a function name or function handle')
  end
  if (~ isa (nboot, 'numeric'))
    error ('bootnhst: NBOOT must be numeric');
  end
  if (any (nboot ~= abs (fix (nboot))))
    error ('bootnhst: NBOOT must contain positive integers')
  end
  if (numel (nboot) > 2)
    error ('bootnhst: A vector for NBOOT cannot have length > 2')
  elseif (numel (nboot) < 2)
    nboot = cat (2, nboot, 0);
  end
  if (nboot(1) < 999)
    error ('bootnhst: The minimum allowable value of NBOOT(1) is 999')
  end 
  if (nboot(2) == 0) && (nvar > 1)
    error (cat (2, 'bootnhst: Jackknife currently only available for', ...
                   ' analysis of univariate data.'))
  end
  if ((nboot(2) == 0) && (~ strcmp (func2str (bootfun), 'mean')))
    if (~ exist ('jackknife','file'))
      if (ISOCTAVE); 
        warning ('bootnhst:jackfail', cat (2, '''jackknife'' function from', ...
                 ' statistics package not found. nboot(2) set to 100.'))
      else
        warning ('bootnhst:jackfail',cat (2, '''jackknife'' function from', ...
                 ' statistics toolbox not found. nboot(2) set to 99.'))
      end
      nboot(2) = 99;
    end
  end
  
  % Error checking
  if (~ isempty (ref) && strcmpi (ref,'pairwise'))
    ref = [];
  end
  if (nargout > 3)
    error ('bootnhst: Only supports up to 3 output arguments')
  end
  if (~ islogical (DisplayOpt) || (numel (DisplayOpt) > 1))
    error ('bootnhst: The value DISPLAYOPT must be a logical scalar value')
  end

  % Data or group exclusion using NaN 
  if (isnumeric (group))
    if (any (isnan (group)))
      data(isnan (group),:) = [];
      group(isnan (group)) = [];
    end
  end
  if (any (any (isnan (data), 2)))
    group(any (isnan (data), 2)) = [];
    data(any (isnan (data), 2),:) = [];
  end

  % Assign non-zero numbers to group labels
  [gnames,jnk,g] = unique (group);
  clear jnk;
  gk = unique (g);
  k = numel (gk);
  if (k > 1)
    if (~ isempty (ref))
      if (isnumeric (ref))
        ref = gk(ismember (gnames, ref));
      else
        ref = gk(strcmp (gnames, ref));
      end
    end
  else
    error ('bootnhst: GROUP must define atleast two groups')
  end
  N = numel (g);
  
  % If applicable, check we have parallel computing capabilities
  if (paropt.UseParallel)
    if (ISOCTAVE)
      software = pkg ('list');
      names = cellfun (@(S) S.name, software, 'UniformOutput', false);
      status = cellfun (@(S) S.loaded, software, 'UniformOutput', false);
      index = find (~ cellfun (@isempty, regexpi (names, '^parallel')));
      if (~ isempty (index))
        if (logical (status{index}))
          PARALLEL = true;
        else
          PARALLEL = false;
        end
      else
        PARALLEL = false;
      end
    else
      if (ismember ('Parallel Computing Toolbox', {info.Name}))
        PARALLEL = true;
      else
        PARALLEL = false;
      end
    end
  end
  
  % If applicable, setup a parallel pool (required for MATLAB)
  if (~ ISOCTAVE)
    % MATLAB
    if (paropt.UseParallel)
      % PARALLEL
      if (paropt.nproc > 0) 
        % MANUAL
        try 
          pool = gcp ('nocreate'); 
          if (isempty (pool))
            if (paropt.nproc > 1)
              % Start parallel pool with nproc workers
              parpool (paropt.nproc);
            else
              % Parallel pool is not running and nproc is 1 so run function
              % evaluations in serial
              paropt.UseParallel = false;
            end
          else
            if (pool.NumWorkers ~= paropt.nproc)
              % Check if number of workers matches nproc and correct it
              % accordingly if not
              delete (pool);
              if (paropt.nproc > 1)
                parpool (paropt.nproc);
              end
            end
          end
        catch
          % MATLAB Parallel Computing Toolbox is not installed
          warning ('bootnhst:parallel', cat (2, 'Parallel Computing', ...
               ' Toolbox is not installed. Falling back to serial processing.'))
          paropt.UseParallel = false;
          paropt.nproc = 1;
        end
      else
        % AUTOMATIC
        try 
          pool = gcp ('nocreate'); 
          if (isempty (pool))
            % Parallel pool not running, start parallel pool using all
            % available workers
            parpool;
          else
            % Parallel pool is already running, set nproc to the
            % number of workers
            paropt.nproc = pool.NumWorkers;
          end
        catch
          % Parallel toolbox not installed, run function evaluations in serial
          paropt.UseParallel = false;
        end
      end
    end
  else
    if (paropt.UseParallel && (paropt.nproc > 1) && ~PARALLEL)
      if (ISOCTAVE)
        % OCTAVE Parallel Computing Package is not installed or loaded
        warning ('bootnhst:parallel', cat (2, 'Parallel Computing Package', ...
         ' is not installed and/or loaded. Falling back to serial processing.'))
      else
        % MATLAB Parallel Computing Toolbox is not installed or loaded
        warning ('bootnhst:parallel', cat (2, 'Parallel Computing Toolbox', ...
         ' is not installed and/or loaded. Falling back to serial processing.'))
      end
      paropt.UseParallel = false;
      paropt.nproc = 0;
    end
  end

  % Create maxstat anonymous function for bootstrap
  func = @(data) localfunc.maxstat (data, g, nboot(2), bootfun, ref, ISOCTAVE);

  % Perform resampling and calculate bootstrap statistics to estimate sampling
  % distribution under the null hypothesis
  boot (1, 1, false, 1); % set random seed to make resampling deterministic
  % Use newer, faster and balanced (less biased) resampling
  % functions (boot and bootknife)
  if (paropt.UseParallel)
    [jnk, Q] = bootknife (data, nboot(1), func, ...
                           NaN, [], paropt.nproc, [], [], ISOCTAVE);
  else
    [jnk, Q] = bootknife (data, nboot(1), func, NaN, [], 0, [], [], ISOCTAVE);
  end
  
  % Compute the estimate (theta) and it's pooled (weighted mean) sampling
  % variance 
  theta = zeros (k, 1);
  SE = zeros (k, 1);
  Var = zeros (k, 1);
  nk = zeros (size(gk));
  for j = 1:k
    if (nboot(2) == 0)
      nk(j) = sum (g == gk(j));
      if (strcmp (func2str (bootfun), 'mean'))
        theta(j) = mean (data(g == gk(j), :));
        % Quick analytical calculation for the standard error of the mean
        SE(j) = std (data(g == gk(j), :), 0) / sqrt (nk(j));
        if (j == 1); se_method = 'Calculated without resampling'; end;
      else
        theta(j) = bootfun (data(g == gk(j), :));
        % If requested, compute unbiased estimates of the standard error using
        % jackknife resampling
        jackstat = jackknife (bootfun, data(g == gk(j), :));
        SE(j) = sqrt ((nk(j) - 1) / nk(j) * ...
                sum (((mean (jackstat) - jackstat)).^2));
        if (j == 1); se_method = 'Leave-one-out jackknife'; end;
      end
    else
      % Compute unbiased estimate of the standard error by balanced bootknife
      % resampling. Bootknife resampling involves less computation than
      % Jackknife when sample sizes get larger
      theta(j) = bootfun (data(g == gk(j), :));
      nk(j) = sum (g == gk(j));
      bootout = bootknife (data(g == gk(j), :), [nboot(2), 0], bootfun, ...
                           NaN, [], 0, [], [], ISOCTAVE, false);
      SE(j) = bootout.std_error;
      if (j==1); se_method = 'Balanced, bootknife resampling'; end;
    end
    Var(j) = ((nk(j) - 1) / (N - k)) * SE(j)^2;
  end
  if (any (SE == 0))
    error ('bootnhst: Samples must have non-zero standard error')
  end
  if (any (isnan (SE)))
    error (cat (2, 'bootnhst: Evaluating bootfun on the bootknife', ...
                   ' resamples created NaN values for the standard error'))
  end
  nk_bar = sum (nk.^2) ./ sum (nk); % weighted mean sample size
  Var = sum (Var .* nk / nk_bar);   % weighted pooled sampling variance
  df = sum (nk) - k;                % degrees of freedom (if bootfun is @mean)

  % Calculate weights to correct for unequal sample size  
  % when calculating standard error of the difference
  w = nk_bar ./ nk;

  % Prepare to make symmetrical bootstrap-t confidence intervals and
  % 2-tailed p-values, Create empirical distribution function
  [Q, F, P] = bootcdf (Q, true, 1);

  % Compute resolution limit of the p-values as determined by resampling
  % with nboot(1) resamples
  res_lim = 1 / (nboot(1) + 1);

  % Calculate p-values for comparisons adjusted to simultaneously
  % control the Family-Wise Error Rate (FWER)
  if (isempty (ref))
    % Single-step maxT procedure for pairwise comparisons is a resampling
    % version of Tukey-Kramer Honest Significant Difference (HSD) test
    A = ones (k, 1) * gk';
    B = tril (gk * ones (1, k), -1);
    M = [A(:) B(:)];
    ridx = (M(:,2) == 0);
    M(ridx, :) = [];
    n = size (M, 1);
    c = zeros (n, 9);
    c(:,1:2) = M;
    for i = 1:n
      c(i,3) = theta(c(i,1));
      c(i,4) = theta(c(i,2));
      c(i,5) = c(i,4) - c(i,3);
      SED = sqrt (Var * (w(c(i,1)) + w(c(i,2))));
      c(i,6) = abs (c(i,5)) / SED;
      if (c(i,6) < Q(1))
        c(i,7) = interp1 (Q, P, c(i,6), 'linear', 1);
      else
        c(i,7) = interp1 (Q, P, c(i,6), 'linear', res_lim);
      end
      c(i,8) = c(i,5) - SED * interp1 (F, Q, 1 - alpha, 'linear', max (Q));
      c(i,9) = c(i,5) + SED * interp1 (F, Q, 1 - alpha, 'linear', max (Q));
    end
  else
    % Single-step maxT procedure for treatment vs control comparisons is
    % a resampling version of Dunnett's test
    c = zeros (k, 9);
    c(:,1) = ref;
    c(:,3) = theta (ref);
    for j = 1:k
      c(j,2) = gk(j);
      c(j,4) = theta(c(j,2));
      c(j,5) = c(j,4) - c(j,3); 
      SED = sqrt (Var * (w(c(j,1)) + w(c(j,2))));
      c(j,6) = abs (c(j,5)) / SED;
      if (c(j,6) < Q(1))
        c(j,7) = interp1 (Q, P, c(j,6), 'linear', 1);
      else
        c(j,7) = interp1 (Q, P, c(j,6), 'linear', res_lim);
      end
      c(j,8) = c(j,5) - SED * interp1 (F, Q, 1 - alpha, 'linear', max (Q));
      c(j,9) = c(j,5) + SED * interp1 (F, Q, 1 - alpha, 'linear', max (Q));
    end
    c(ref,:) = [];
  end

  % Calculate the (maximum) test statistic and (minimum) p-value for the
  % omnibus test
  maxT = max (c(:,6));
  pval = min (c(:,7));

  % Prepare stats output structure
  stats = struct;
  stats.gnames = gnames;
  stats.ref = ref;
  stats.groups = zeros (k,6);
  stats.groups = zeros (k,6);
  stats.groups(:,1) = gk;
  stats.groups(:,2) = theta;
  stats.groups(:,3) = nk;
  stats.groups(:,4) = SE;
  stats.groups(:,5) = theta - sqrt ((0.5 * (w + 1)) .* Var / 2) * ...
                              interp1 (F, Q, 1 - alpha, 'linear', max (Q));
  stats.groups(:,6) = theta + sqrt ((0.5 * (w + 1)) .* Var / 2) * ...
                              interp1 (F, Q, 1 - alpha, 'linear', max (Q));
  stats.Var = Var;
  stats.maxT = maxT;
  if (any (strcmp (func2str (bootfun), {'mean','smoothmedian'})))
    stats.df = df;
  else
    stats.df = [];
  end
  stats.nboot = nboot;
  stats.alpha = alpha;
  stats.bootstat = Q;

  % Print output and plot graph with confidence intervals if no output
  % arguments are requested
  cols = [1,2,5,6,7]; % columns in c that we want to print data for
  if ((nargout == 0) || (DisplayOpt == true))
    if (~ iscellstr (gnames))
      gnames = cellstr (num2str (gnames));
    end
    fprintf (cat (2, '\nSummary of bootstrap null hypothesis (H0)', ...
                     ' significance test(s)\n', ...
                     '*******************************************', ...
                     '***********************************\n'));
    fprintf ('Bootstrap settings: \n');
    fprintf (' Function: %s\n', bootfun_str);
    fprintf (' Bootstrap resampling method: Balanced, bootknife resampling\n')
    fprintf (' Number of bootstrap resamples: %u \n', nboot(1));
    fprintf (' Method for estimating standard errors: %s\n', se_method)
    if (nboot(2) > 0)
      fprintf (cat (2, ' Number of bootknife resamples used to estimate', ...
                       ' standard errors: %u \n'), nboot(2));
    end
    if (isempty (ref))
      fprintf (' Multiple comparison method:%s \n', ... 
               ' Single-step maxT procedure based on Tukey-Kramer');
    else
      fprintf (' Multiple comparison method:%s \n', ...
               ' Single-step maxT procedure based on Dunnett');
      fprintf (' Reference group used for comparisons: %s \n', gnames{ref});
    end
    fprintf (cat (2, '------------------------------------------------', ...
                     '------------------------------\n\n'));
    if (isempty (ref))
      fprintf (cat (2, 'Overall hypothesis test from single-step maxT', ...
                       ' procedure\nH0: Groups of data are all sampled', ... 
                       ' from the same population\n\n'));
    else
      fprintf (cat (2, 'Overall hypothesis test from single-step maxT', ... 
                       ' procedure\nH0: Groups of data are all sampled', ...
                       ' from the same population as data in REF\n\n'));
    end
    if (any (strcmp (func2str (bootfun), {'mean','smoothmedian'})))
      dfstr = sprintf ('Degrees of freedom = %u\n', df);
    else
      dfstr = '';
    end
    if (pval <= 0.001)
      fprintf (cat (2, 'Maximum t = %.2f, p = <.001 \n', dfstr, ...
                '-----------------------------------------------', ...
                '-------------------------------\n'), maxT);
    elseif (pval > 0.999)
      fprintf (cat (2, 'Maximum t = %.2f, p = 1.000 \n', dfstr, ...
                '-----------------------------------------------', ...
                '-------------------------------\n'), maxT);
    else
      fprintf (cat (2, 'Maximum t = %.2f, p = .%03u \n', dfstr, ...
                '-----------------------------------------------', ...
                '-------------------------------\n'), maxT, round (pval*1000));
    end
    if (size (c,1) >= 1)
      fprintf (cat (2, 'POST HOC TESTS with control of the FWER by the', ...
                       ' single-step maxT procedure\n', ...
                       '----------------------------------------------', ...
                       '--------------------------------\n', ...
                       '| Comparison |  Reference # |       Test # |', ...
                       '  Difference |       t |       p |\n', ...
                       '|------------|--------------|--------------|', ...
                       '-------------|---------|---------|\n'), df);
      if (isempty (ref))
        for i = 1:n
          tmp = num2cell (c(i, cols));
          tmp{end} = round (tmp{end} * 1000);
          if (c(i,7) <= 0.001)
            tmp(end) = [];
            fprintf (cat (2, '| %10u | %12u | %12u | %+11.2e | %7.2f |', ...
                             '   <.001 |***\n'), i, tmp{:});
          elseif (c(i,7) > 0.999)
            tmp(end) = [];
            fprintf (cat (2, '| %10u | %12u | %12u | %+11.2e | %7.2f |', ...
                             '   1.000 |\n'), i, tmp{:});
          else
            fprintf (cat (2, '| %10u | %12u | %12u | %+11.2e | %7.2f |', ...
                             '    .%03u |'), i, tmp{:});
            if (c(i,7) < 0.01)
              fprintf ('**\n')
            elseif (c(i,7) < 0.05)
              fprintf ('*\n')
            else
              fprintf ('\n')
            end
          end
        end
      else
        for j = 1:k-1
          tmp = num2cell (c(j, cols));
          tmp{end} = round (tmp{end} * 1000);
          if (c(j,7) <= 0.001)
            tmp(end) = [];
            fprintf (cat (2, '| %10u | %12u | %12u | %+11.2e | %7.2f |', ...
                             '   <.001 |***\n'), j, tmp{:});
          elseif (c(j,7) > 0.999)
            tmp(end) = [];
            fprintf (cat (2, '| %10u | %12u | %12u | %+11.2e | %7.2f |', ...
                             '   1.000 |\n'), j, tmp{:});
          else
            fprintf (cat (2, '| %10u | %12u | %12u | %+11.2e | %7.2f |', ...
                             '    .%03u |'), j, tmp{:});
            if (c(j,7) < 0.01)
              fprintf ('**\n')
            elseif (c(j,7) < 0.05)
              fprintf ('*\n')
            else
              fprintf ('\n')
            end
          end
        end
      end
      fprintf (cat (2, '\n-----------------------------------------------', ...
                      '-------------------------------\n', ...
                      '|    GROUP # |                                   ', ...
                      '                GROUP label |\n', ...
                      '|------------|-----------------------------------', ...
                      '----------------------------|\n'));
      for j = 1:k
        fprintf ('| %10u | %61s |\n', gk(j), gnames{j});
      end
      fprintf ('\n')
    end

    % Plot graph of the difference in bootfun for each comparison with
    % 100*(1-alpha)% confidence intervals
    figure;
    nc = size(c,1);                    % Calculate number of comparisons to plot
    plot ([0; 0], [0; nc + 1]', 'k:'); % Plot vertical dashed line at 0 effect
    set (gca, 'Ydir', 'reverse')       % Flip y-axis direction
    ylim ([0.5, nc + 0.5]);            % Set y-axis limits
    hold on
    for i = 1 : nc
      if (c(i,7) < 0.05)
        % Plot marker for the difference estimate
        plot (c(i, 5), i, 'or', 'MarkerFaceColor', 'r');
        % Plot line for each confidence interval
        plot ([c(i, 8), c(i, 9)], i * ones (2, 1), 'r-');
      else
        % Plot marker for the difference estimate
        plot (c(i,5), i, 'ob', 'MarkerFaceColor', 'b');
        % Plot line for each confidence interval 
        plot ([c(i,8), c(i,9)], i * ones (2, 1), 'b-');   
      end
    end
    hold off
    xlabel (sprintf (cat (2, '%g%% bootstrap-t confidence interval for the', ...
                             ' difference'), 100 * (1 - alpha)));
    ylabel ('Comparison number (Test - Reference)');   

  end

end

%--------------------------------------------------------------------------

function maxT = maxstat (Y, g, nboot, bootfun, ref, ISOCTAVE)

  % Helper function file required for bootnhst
  % Calculate maximum test statistic
  
  % maxstat cannot be a subfunction or nested function since 
  % Octave parallel threads won't be able to find it

  % Calculate the size of the data (N) and the number (k) of unique groups
  N = size (g, 1);
  gk = unique (g);
  k = numel (gk);

  % Compute the estimate (theta) and it's pooled (weighted mean) sampling
  % variance 
  theta = zeros (k, 1);
  SE = zeros (k, 1);
  Var = zeros (k, 1);
  nk = zeros (size (gk));
  for j = 1:k
    if (nboot == 0)
      nk(j) = sum (g == gk(j));
      if strcmp (func2str (bootfun), 'mean')
        theta(j) = mean (Y(g == gk(j), :));
        % Quick calculation for the standard error of the mean
        SE(j) = std (Y(g == gk(j), :), 0) / sqrt (nk(j));
      else
        theta(j) = bootfun (Y(g == gk(j), :));
        % If requested, compute unbiased estimates of the standard error
        % using jackknife resampling
        jackstat = jackknife (bootfun, Y(g == gk(j), :));
        SE(j) = sqrt ((nk(j) - 1) / nk(j) ...
                * sum (((mean (jackstat) - jackstat)).^2));
      end
    else
      % Compute unbiased estimate of the standard error by balanced bootknife
      % resampling. Bootknife resampling involves less computation than
      % Jackknife when sample sizes get larger
      theta(j) = bootfun (Y(g == gk(j), :));
      nk(j) = sum (g == gk(j));
      bootout = bootknife (Y(g == gk(j), :), [nboot, 0], bootfun, ...
                           NaN, [], 0, [], [], ISOCTAVE, false);
      SE(j) = bootout.std_error;
    end
    Var(j) = ((nk(j) - 1) / (N - k)) * SE(j)^2;
  end
  if (any (isnan (SE)))
    error (cat (2, 'bootnhst:maxstat: Evaluating bootfun on the bootknife', ...
                   ' resamples created NaN values for the standard error'))
  end
  nk_bar = sum (nk.^2) ./ sum (nk);  % weighted mean sample size
  Var = sum (Var .* nk / nk_bar);    % weighted pooled sampling variance

  % Calculate weights to correct for unequal sample size  
  % when calculating standard error of the difference
  w = nk_bar ./ nk;

  % Calculate the maximum test statistic 
  if (isempty (ref))
    % Calculate Tukey-Kramer test statistic (without sqrt(2) factor)
    %
    % Bibliography:
    %  [1] en.wikipedia.org/wiki/Tukey%27s_range_test
    %  [2] cdn.graphpad.com/faq/1688/file/MulitpleComparisonAlgorithmsPrism8.pdf
    %  [3] www.graphpad.com/guides/prism/
    %          latest/statistics/stat_the_methods_of_tukey_and_dunne.htm
    idx = logical (triu (ones (k, k), 1));
    i = (1 : k)' * ones (1, k);
    j = ones (k, 1) * (1 : k);
    t = abs (theta(i(idx)) - theta(j(idx))) ./ ...
        sqrt (Var * (w(i(idx)) + w(j(idx))));
  else
    % Calculate Dunnett's test statistic 
    t = abs ((theta - theta(ref))) ./ sqrt (Var * (w + w(ref)));
  end
  maxT = max(t);
  
end

%--------------------------------------------------------------------------

%!demo
%!
%! ## EXAMPLE 1A: 
%! ## ONE-WAY ANOVA WITH EQUAL SAMPLE SIZES: Treatment vs. Control (1)
%!
%! y = [111.39 110.21  89.21  76.64  95.35  90.97  62.78;
%!      112.93  60.36  92.29  59.54  98.93  97.03  79.65;
%!       85.24 109.63  64.93  75.69  95.28  57.41  75.83;
%!      111.96 103.40  75.49  76.69  77.95  93.32  78.70];
%! g = [1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7];
%!
%! bootnhst (y(:),g(:),'ref',1,'nboot',4999);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## EXAMPLE 1B: 
%! ## ROBUST ONE-WAY ANOVA WITH EQUAL SAMPLE SIZES: Treatment vs. Control (1)
%!
%! y = [111.39 110.21  89.21  76.64  95.35  90.97  62.78;
%!      112.93  60.36  92.29  59.54  98.93  97.03  79.65;
%!       85.24 109.63  64.93  75.69  95.28  57.41  75.83;
%!      111.96 103.40  75.49  76.69  77.95  93.32  78.70];
%! g = [1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7];
%!
%! bootnhst (y(:), g(:), 'ref', 1, 'nboot', 4999, 'bootfun', 'robust');
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## EXAMPLE 2A:
%! ## COMPARISON OF TWO INDEPENDENT GROUPS WITH UNEQUAL SAMPLE SIZES 
%! ## (analagous to Student's t-test)
%!
%! y =    [54       43
%!         23       34
%!         45       65
%!         54       77
%!         45       46
%!        NaN       65];
%! g = {'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'};
%!
%! bootnhst (y(:), g(:), 'ref', 'male', 'nboot', 4999);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## EXAMPLE 2B:
%! ## ONE-WAY ANOVA WITH UNEQUAL SAMPLE SIZES: pairwise comparisons
%!
%! y = [54  87  45
%!      23  98  39
%!      45  64  51
%!      54  77  49
%!      45  89  50
%!      47 NaN  55];
%! g = [ 1   2   3
%!       1   2   3
%!       1   2   3 
%!       1   2   3
%!       1   2   3
%!       1   2   3];
%!
%! bootnhst (y(:), g(:), 'nboot', 4999);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## EXAMPLE 2C:
%! ## COMPARE STANDARD DEVIATIONS BETWEEN 3 GROUPS: pairwise comparisons
%!
%! y = [54  87  45
%!      23  98  39
%!      45  64  51
%!      54  77  49
%!      45  89  50
%!      47 NaN  55];
%! g = [ 1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3];
%! p = bootnhst (y(:),g(:),'bootfun',{@std,1});

%!demo
%!
%! ## EXAMPLE 3:
%! ## COMPARE CORRELATION COEFFICIENTS BETWEEN 2 DATA SETS
%! Y = randn (20, 2); g = [zeros(10, 1); ones(10, 1)];
%! func = @(M) cor (M(:,1), M(:,2));
%!
%! bootnhst (Y, g, 'bootfun', func);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## EXAMPLE 4:
%! ## COMPARE SLOPES FROM LINEAR REGRESSION ON 2 DATA SETS
%! y = randn (20, 1); x = randn (20, 1); X = [ones(20, 1), x];
%! g = [zeros(10, 1); ones(10, 1)];
%! func = @(M) subsref (M(:,2:end) \ M(:,1), ...
%!                      struct ('type', '()', 'subs', {{2}}));
%!
%! bootnhst ([y, X], g, 'bootfun', func);
%!
%! ## Please be patient, the calculations will be completed soon...

%!test
%! y = [111.39 110.21  89.21  76.64  95.35  90.97  62.78;
%!      112.93  60.36  92.29  59.54  98.93  97.03  79.65;
%!       85.24 109.63  64.93  75.69  95.28  57.41  75.83;
%!      111.96 103.40  75.49  76.69  77.95  93.32  78.70];
%! g = [1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7];
%! p = bootnhst (y(:),g(:),'ref',1,'nboot',[999,0],'DisplayOpt',false);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   # test boot m-file result
%!   assert (p, 0.01264993609963437, 1e-06);
%! end
%! p = bootnhst (y(:),g(:),'nboot',[999,0],'DisplayOpt',false);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   # test boot m-file result
%!   assert (p, 0.04197298481251393, 1e-06);
%! end
%! # Result from anova1 is 0.0387

%!test
%! y = [54       43
%!      23       34
%!      45       65
%!      54       77
%!      45       46
%!     NaN       65];
%! g = {'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'};
%! p = bootnhst (y(:),g(:),'ref','male','nboot',[999,0],'DisplayOpt',false);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   # test boot m-file result
%!   assert (p, 0.2775241109820157, 1e-06);
%! end
%! p = bootnhst (y(:),g(:),'nboot',[999,0],'DisplayOpt',false);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   # test boot m-file result
%!   assert (p, 0.2775241109820157, 1e-06);
%! end
%! # Result from anova1 is 0.2613

%!test
%! y = [54  87  45
%!      23  98  39
%!      45  64  51
%!      54  77  49
%!      45  89  50
%!      47 NaN  55];
%! g = [ 1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3];
%! p = bootnhst (y(:),g(:),'nboot',[999,0],'DisplayOpt',false);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   # test boot m-file result
%!   assert (p, 0.001, 1e-06); # truncated at 0.001
%! end
%! # Result from anova1 is 4.162704768129188e-05

%!test
%! y = [54  87  45
%!      23  98  39
%!      45  64  51
%!      54  77  49
%!      45  89  50
%!      47 NaN  55];
%! g = [ 1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3];
%! p = bootnhst (y(:),g(:),'bootfun',@(y)std(y,1),'DisplayOpt',false);
%! p = bootnhst (y(:),g(:),'bootfun',{@std,1},'DisplayOpt',false);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   # test boot m-file result
%!   assert (p, 0.4421449140843043, 1e-06);
%! end

%!test
%! % Compare correlation coefficients
%! Y = randn (20, 2); g = [zeros(10, 1); ones(10, 1)];
%! func = @(M) cor (M(:,1), M(:,2));
%! p = bootnhst (Y, g, 'bootfun', func, 'DisplayOpt', false);

%!test
%! % Compare slopes from linear regression
%! y = randn (20, 1); x = randn (20, 1); X = [ones(20, 1), x];
%! g = [zeros(10, 1); ones(10, 1)];
%! func = @(M) subsref (M(:,2:end) \ M(:,1), ...
%!                      struct ('type', '()', 'subs', {{2}}));
%! p = bootnhst ([y, X], g, 'bootfun', func, 'DisplayOpt', false);