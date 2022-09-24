%  Function File: bootknife
%
%  Bootknife (bootstrap) resampling
%
%  This function takes a data sample (containing n rows) and uses bootstrap 
%  techniques to calculate a bias of the parameter estimate, a standard 
%  error, and 95% confidence intervals [1]. Specifically, the method uses 
%  bootknife resampling [2], which involves creating leave-one-out jackknife
%  samples of size n - 1 and then drawing samples of size n with replacement
%  from the jackknife samples. The resampling of data rows is balanced in order
%  to reduce Monte Carlo error [3,4]. The confidence intervals that are returned
%  are percentile bootstrap confidence intervals.
%
%  stats = bootknife(data)
%  stats = bootknife(data,nboot)
%  stats = bootknife(data,nboot,bootfun)
%  stats = bootknife(data,nboot,{bootfun,bootfun_args})
%  stats = bootknife(data,nboot,bootfun,alpha)
%  stats = bootknife(data,nboot,bootfun,alpha,strata)
%  stats = bootknife(data,nboot,bootfun,alpha,strata,nproc)
%  stats = bootknife(data,2000,@mean,0.05,[],0)      % Default values
%  [stats,bootstat] = bootknife(...)
%  [stats,bootstat] = bootknife(...)
%  [stats,bootstat,bootsam] = bootknife(...)
%  bootknife(data,...);
%
%  stats = bootknife(data) resamples from the rows of a data sample (column
%  vector or a matrix) and returns a structure with the following fields:
%    original: contains the result of applying bootfun to the data (x)
%    bias: contains the bootstrap estimate of bias
%    std_error: contains the bootstrap standard error
%    CI_lower: contains the lower bound of the bootstrap confidence interval
%    CI_upper: contains the upper bound of the bootstrap confidence interval
%  By default, the statistics relate to bootfun being @mean and the confidence
%  intervals are 95% percentile bootstrap intervals. 
%
%  stats = bootknife(data,nboot) also specifies the number of bootstrap
%  samples. nboot must be scalar, positive integer. By default, nboot is 2000.
%
%  stats = bootknife(data,nboot,bootfun) specifies bootfun, a function handle,
%  a string indicating the name of the function to apply to the data (and each
%  bootstrap resample), or a cell array where the first cell is the function
%  handle or string, and other cells being arguments for that function, where
%  the function must take data for the first input argument. bootfun can return
%  a scalar value or vector. The default value(s) of bootfun is/are the (column)
%  mean(s).
%    Note that bootfun MUST calculate a statistic representative of the
%  finite data sample, it should NOT be an estimate of a population
%  parameter. For example, for the variance, set bootfun to {@var,1}, not
%  @var or {@var,0}. Smooth functions of the data are preferable, (e.g. use
%  smoothmedian function instead of ordinary median).
%
%  stats = bootknife(data,nboot,bootfun,alpha) where alpha sets the lower
%  and upper confidence interval ends to be 100 * (alpha/2)% and 100 *
%  (1-alpha/2)% respectively. Central coverage of the intervals is thus
%  100*(1-alpha)%. alpha should be a scalar value between 0 and 1. If alpha
%  is empty, NaN is returned for the confidence interval ends. The default
%  alpha is 0.05. 
%
%  stats = bootknife(data,nboot,bootfun,alpha,strata) also sets strata,
%  which are identifiers that define the grouping of the data rows
%  for stratified bootstrap resampling. strata should be a column vector
%  or cell array the same number of rows as the data. When resampling is
%  stratified, the groups (or stata) of data are equally represented across
%  the bootstrap resamples. If this input argument is not specified or is
%  empty, no stratification of resampling is performed. 
%
%  stats = bootknife(data,nboot,bootfun,alpha,strata,nproc) sets the
%  number of parallel processes to use to accelerate computations on
%  multicore machines, specifically non-vectorized function evaluations,
%  double bootstrap resampling and jackknife function evaluations. This
%  feature requires the Parallel package (in Octave), or the Parallel
%  Computing Toolbox (in Matlab).
%
%  [stats,bootstat] = bootknife(...) also returns bootstat, a vector of
%  statistics calculated over the (first, or outer layer of) bootstrap
%  resamples. 
%
%  [stats,bootstat,bootsam] = bootknife(...) also returns bootsam, the
%  matrix of indices (32-bit integers) used for the (first, or outer
%  layer of) bootstrap resampling. Each column in bootsam corresponds
%  to one bootstrap resample and contains the row indices of the values
%  drawn from the nonscalar data argument to create that sample.
%
%  bootknife(data,...); returns a pretty table of the output including
%  the bootstrap settings and the result of evaluating bootfun on the
%  data along with bootstrap estimates of bias, standard error, and
%  lower and upper 100*(1-alpha)% confidence limits.
%
%  Requirements: The function file boot.m (or better boot.mex) also
%  distributed in the statistics-bootstrap package. The 'robust' option
%  for bootfun requires smoothmedian.m (or better smoothmedian.mex).
%
%  Bibliography:
%  [1] Efron, and Tibshirani (1993) An Introduction to the
%        Bootstrap. New York, NY: Chapman & Hall
%  [2] Hesterberg T.C. (2004) Unbiasing the Bootstrap???Bootknife Sampling
%        vs. Smoothing; Proceedings of the Section on Statistics & the
%        Environment. Alexandria, VA: American Statistical Association.
%  [3] Davison et al. (1986) Efficient Bootstrap Simulation.
%        Biometrika, 73: 555-66
%  [4] Gleason, J.R. (1988) Algorithms for Balanced Bootstrap Simulations.
%        The American Statistician. Vol. 42, No. 4 pp. 263-266
%
%  bootknife v1.8.0.0 (23/09/2022)
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


function [stats, bootstat, BOOTSAM] = bootknife (x, nboot, bootfun, alpha, strata, ncpus, REF, ISOCTAVE, BOOTSAM)
  
  % Error checking
  if (nargin < 1)
    error ('data must be provided')
  end
  if ~(size(x, 1) > 1)
    error ('data must contain more than one row')
  end
  
  % Set defaults or check for errors
  if (nargin < 2)
    nboot = 2000;
  else
    if isempty(nboot)
      nboot = 2000;
    else
      if ~isa (nboot, 'numeric')
        error('nboot must be numeric');
      end
      if any (nboot ~= abs (fix (nboot)))
        error ('nboot must contain a positive integer')
      end    
      if (numel (nboot) > 1)
        error ('nboot must be a scalar value')
      end
    end
  end
  if (nargin < 3)
    bootfun = @mean;
  else
    if iscell(bootfun)
      func = bootfun{1};
      args = bootfun(2:end);
      bootfun = @(x) func (x, args{:});
    end
    if ischar (bootfun)
      if strcmpi(bootfun,'robust')
        bootfun = 'smoothmedian';
      end
      % Convert character string of a function name to a function handle
      bootfun = str2func (bootfun);
    end
    if ~isa (bootfun, 'function_handle')
      error('bootfun must be a function name or function handle');
    end
  end
  if (nargin < 4)
    alpha = 0.05;
  elseif ~isempty (alpha) 
    if ~isa (alpha,'numeric') || numel (alpha)~=1
      error('alpha must be a numeric scalar value');
    end
    if (alpha <= 0) || (alpha >= 1)
      error('alpha must be a value between 0 and 1');
    end
  end
  if (nargin < 5)
    strata = [];
  elseif ~isempty (strata) 
    if size (strata, 1) ~= size (x, 1)
      error('strata should be a column vector or cell array with the same number of rows as the data')
    end
  end
  if (nargin < 6)
    ncpus = 0;    % Ignore parallel processing features
  elseif ~isempty (ncpus) 
    if ~isa (ncpus, 'numeric')
      error('ncpus must be numeric');
    end
    if any (ncpus ~= abs (fix (ncpus)))
      error ('ncpus must be a positive integer')
    end    
    if (numel (ncpus) > 1)
      error ('ncpus must be a scalar value')
    end
  end
  % REF, ISOCTAVE and BOOTSAM are undocumented input arguments required for some of the functionality of bootknife
  if (nargin < 8)
    % Check if running in Octave (else assume Matlab)
    info = ver; 
    ISOCTAVE = any (ismember ({info.Name}, 'Octave'));
  end
  if ISOCTAVE
    ncpus = min(ncpus, nproc);
  else
    ncpus = min(ncpus, feature('numcores'));
  end

  % Determine properties of the data (x)
  [n, nvar] = size (x);

  % Evaluate bootfun on the data
  T0 = bootfun (x);
  if all (size (T0) > 1)
    error ('bootfun must return either a scalar or a vector')
  end

  % If data is univariate, check whether bootfun is vectorized
  if (nvar == 1)
      try
        chk = bootfun (cat (2,x,x));
        if ( all (size (chk) == [1, 2]) && all (chk == bootfun (x)) )
          vectorized = true;
        else
          vectorized = false;
        end
      catch
        vectorized = false;
      end
  else
    vectorized = false;
  end
  
  % If applicable, check we have parallel computing capabilities
  if (ncpus > 1)
    if ISOCTAVE  
      pat = '^parallel';
      software = pkg('list');
      names = cellfun(@(S) S.name, software, 'UniformOutput', false);
      status = cellfun(@(S) S.loaded, software, 'UniformOutput', false);
      index = find(~cellfun(@isempty,regexpi(names,pat)));
      if ~isempty(index)
        if logical(status{index})
          PARALLEL = true;
        else
          PARALLEL = false;
        end
      else
        PARALLEL = false;
      end
    else
      try
        retval = ~isempty(getCurrentTask()) && (matlabpool('size') > 0);
      catch err
        if ~strcmp(err.identifier, 'MATLAB:UndefinedFunction')
          rethrow(err);
        end
        PARALLEL = false;
      end
    end
  end
  
  % If applicable, setup a parallel pool (required for MATLAB)
  if ~ISOCTAVE
    % MATLAB
    % bootfun is not vectorized
    if (ncpus > 0) 
      % MANUAL
      try 
        pool = gcp ('nocreate'); 
        if isempty (pool)
          if (ncpus > 1)
            % Start parallel pool with ncpus workers
            parpool (ncpus);
          else
            % Parallel pool is not running and ncpus is 1 so run function evaluations in serial
            ncpus = 1;
          end
        else
          if (pool.NumWorkers ~= ncpus)
            % Check if number of workers matches ncpus and correct it accordingly if not
            delete (pool);
            if (ncpus > 1)
              parpool (ncpus);
            end
          end
        end
      catch
        % MATLAB Parallel Computing Toolbox is not installed
        warning('MATLAB Parallel Computing Toolbox is not installed. Falling back to serial processing.')
        ncpus = 1;
      end
    end
  else
    if (ncpus > 1) && ~PARALLEL
      if ISOCTAVE
        % OCTAVE Parallel Computing Package is not installed or loaded
        warning('OCTAVE Parallel Computing Package is not installed and/or loaded. Falling back to serial processing.')
      else
        % MATLAB Parallel Computing Toolbox is not installed or loaded
        warning('MATLAB Parallel Computing Toolbox is not installed and/or loaded. Falling back to serial processing.')
      end
      ncpus = 0;
    end
  end

  % If the function of the data is a vector, calculate the statistics for each element 
  sz = size (T0);
  m = prod (sz);
  if (m > 1)
    if (nvar == m)
      try
        % If data is multivariate, check whether bootfun is vectorized
        % Bootfun will be evaluated for each column of x, considering each of them as univariate data vectors
        chk = bootfun (cat (2,x(:,1),x(:,1)));
        if ( all (size (chk) == [1, 2]) && all (chk == bootfun (x(:,1))) )
          vectorized = true;
        end
      catch
        % Do nothing
      end
    end
    % Use bootknife for each element of the output of bootfun
    % Note that row indices in the resamples are the same for all columns of data
    stats = struct ('original',zeros(sz),...
                    'bias',zeros(sz),...
                    'std_error',zeros(sz),...
                    'CI_lower',zeros(sz),...
                    'CI_upper',zeros(sz));
    bootstat = zeros (m, nboot);
    if vectorized
      for j = 1:m
        if j > 1
          [stats(j), bootstat(j,:)] = bootknife (x(:,j), nboot, bootfun, alpha, strata, ncpus, [], ISOCTAVE, BOOTSAM);
        else
          [stats(j), bootstat(j,:), BOOTSAM] = bootknife (x(:,j), nboot, bootfun, alpha, strata, ncpus, [], ISOCTAVE);
        end
      end
    else
      for j = 1:m
        out = @(t) t(j);
        func = @(x) out (bootfun (x));
        if j > 1
          [stats(j), bootstat(j,:)] = bootknife (x, nboot, func, alpha, strata, ncpus, [], ISOCTAVE, BOOTSAM);
        else
          [stats(j), bootstat(j,:), BOOTSAM] = bootknife (x, nboot, func, alpha, strata, ncpus, [], ISOCTAVE);
        end
      end
    end
    % Print output if no output arguments are requested
    if (nargout == 0) 
      print_output(stats);
    end
    return
  end

  % Evaluate strata input argument
  if ~isempty (strata)
    if ~isnumeric (strata)
      % Convert strata to numeric ID
      [junk1,junk2,strata] = unique (strata,'legacy');
      clear junk1 junk2;
    end
    % Get strata IDs
    gid = unique (strata,'legacy');  % strata ID
    K = numel (gid);        % number of strata
    % Create strata matrix
    g = false (n,K);
    for k = 1:K
      g(:, k) = (strata == gid(k));
    end
    nk = sum(g).';          % strata sample sizes
  else 
    g = ones(n,1);
  end

  % Perform balanced bootknife resampling
  if nargin < 9
    if ~isempty (strata)
      if (nvar > 1) || (nargout > 2)
        % If we need BOOTSAM, can save some memory by making BOOTSAM an int32 datatype
        BOOTSAM = zeros (n, nboot, 'int32'); 
        for k = 1:K
          BOOTSAM(g(:, k),:) = boot (find (g(:, k)), nboot, true);
        end
      else
        % For more efficiency, if we don't need BOOTSAM, we can directly resample values of x
        BOOTSAM = [];
        X = zeros (n, nboot);
        for k = 1:K
          X(g(:, k),:) = boot (x(g(:, k),:), nboot, true);
        end
      end
    else
      if (nvar > 1) || (nargout > 2)
        % If we need BOOTSAM, can save some memory by making BOOTSAM an int32 datatype
        BOOTSAM = zeros (n, nboot, 'int32');
        BOOTSAM(:,:) = boot (n, nboot, true);
      else
        % For more efficiency, if we don't need BOOTSAM, we can directly resample values of x
        BOOTSAM = [];
        X = boot (x, nboot, true);
      end
    end
  end
  if isempty(BOOTSAM)
    if vectorized
      % Vectorized evaluation of bootfun on the data resamples
      bootstat = bootfun (X);
    else
      if (ncpus > 1)
        % Evaluate bootfun on each bootstrap resample in PARALLEL
        if ISOCTAVE
          % OCTAVE
          bootstat = parcellfun (ncpus, bootfun, num2cell (X, 1));
        else
          % MATLAB
          bootstat = zeros (1, nboot);
          parfor b = 1:nboot; bootstat(b) = cellfunc (X(:, b)); end;
        end
      else
        bootstat = cellfun (bootfun, num2cell (X, 1));
      end
    end
  else
    if vectorized
      % Vectorized implementation of data sampling (using BOOTSAM) and evaluation of bootfun on the data resamples 
      % Perform data sampling
      X = x(BOOTSAM);
      % Function evaluation on bootknife sample
      bootstat = bootfun (X);
    else 
      cellfunc = @(BOOTSAM) bootfun (x(BOOTSAM, :));
      if (ncpus > 1)
        % Evaluate bootfun on each bootstrap resample in PARALLEL
        if ISOCTAVE
          % OCTAVE
          bootstat = parcellfun (ncpus, cellfunc, num2cell (BOOTSAM, 1));
        else
          % MATLAB
          bootstat = zeros (1, nboot);
          parfor b = 1:nboot; bootstat(b) = cellfunc (BOOTSAM(:, b)); end;
        end
      else
        % Evaluate bootfun on each bootstrap resample in SERIAL
        cellfunc = @(BOOTSAM) bootfun (x(BOOTSAM, :));
        bootstat = cellfun (cellfunc, num2cell (BOOTSAM, 1));
      end
    end
  end

  % Bootstrap bias estimation
  bias = mean (bootstat) - T0;
  % Bootstrap standard error
  se = std (bootstat, 1);
  if ~isempty(alpha)
    % Percentile confidence intervals 
    % Create distribution functions
    [cdf, t1] = empcdf (bootstat, 1);
    ci = arrayfun ( @(p) interp1 (cdf, t1, p, 'linear'), [alpha / 2, 1 - alpha / 2]);
  else
    ci = nan (1, 2);
  end

  % Prepare stats output argument
  stats = struct;
  stats.original = T0;
  stats.bias = bias;
  stats.std_error = se;
  stats.CI_lower = ci(1);
  stats.CI_upper = ci(2);
  
  % Print output if no output arguments are requested
  if (nargout == 0) 
    print_output(stats);
  end

  %--------------------------------------------------------------------------

  function print_output(stats)

      fprintf (['\nSummary of non-parametric bootstrap estimates of bias and precision\n',...
                '******************************************************************************\n\n']);
      fprintf ('Bootstrap settings: \n');
      fprintf (' Function: %s\n',func2str(bootfun));
      fprintf (' Resampling method: Balanced, bootknife resampling \n')
      fprintf (' Number of resamples: %u \n', nboot);
      fprintf (' Confidence interval type: Percentile \n');
      fprintf (' Confidence interval coverage: %g%% \n\n',100*(1-alpha));
      fprintf ('Bootstrap Statistics: \n');
      fprintf (' original       bias           std_error      CI_lower       CI_upper    \n');
      for i = 1:m
        fprintf (' %#-+12.6g   %#-+12.6g   %#-+12.6g   %#-+12.6g   %#-+12.6g \n',... 
                 [stats(i).original, stats(i).bias, stats(i).std_error, stats(i).CI_lower, stats(i).CI_upper]);
      end
      fprintf ('\n');
      
  end

end

%--------------------------------------------------------------------------

function [F, x] = empcdf (bootstat, c)

  % Subfunction to calculate empirical cumulative distribution function of bootstat
  %
  % Set c to:
  %  1 to have a complete distribution with F ranging from 0 to 1
  %  0 to avoid duplicate values in x
  %
  % Unlike ecdf, empcdf uses a denominator of N+1

  % Check input argument
  if ~isa(bootstat,'numeric')
    error('bootstat must be numeric')
  end
  if all(size(bootstat)>1)
    error('bootstat must be a vector')
  end
  if size(bootstat,2)>1
    bootstat = bootstat.';
  end

  % Create empirical CDF
  x = sort(bootstat);
  N = sum(~isnan(bootstat));
  [x,F] = unique(x,'rows','last','legacy');
  F = F/(N+1);

  % Apply option to complete the CDF
  if c > 0
    x = [x(1);x;x(end)];
    F = [0;F;1];
  end

  % Remove impossible values
  F(isnan(x)) = [];
  x(isnan(x)) = [];
  F(isinf(x)) = [];
  x(isinf(x)) = [];

end