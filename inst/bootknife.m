%  Function File: bootknife
%
%  Bootknife (bootstrap) resampling
%
%  This function takes a data sample (containing n rows) and uses bootstrap 
%  techniques to calculate a bias of the parameter estimate, a standard 
%  error, and 95% confidence intervals. Specifically, the method uses 
%  bootknife resampling [1], which involves creating leave-one-out 
%  jackknife samples of size n - 1 and then drawing samples of size n with 
%  replacement from the jackknife samples. The resampling of data rows is 
%  balanced in order to reduce Monte Carlo error [2,3]. By default, the 
%  bootstrap confidence intervals are bias-corrected and accelerated (BCa) 
%  [4-5]. BCa intervals are fast to compute and have good coverage and 
%  correctness when combined with bootknife resampling as it is here [1], 
%  but it may not have the intended coverage when sample size gets very 
%  small. If double bootstrap is requested, the algorithm uses calibration 
%  to improve the accuracy of the bias estimate and confidence intervals
%  for small-to-medium sample sizes [6-8]. 
%
%  stats = bootknife(data)
%  stats = bootknife(data,nboot)
%  stats = bootknife(data,nboot,bootfun)
%  stats = bootknife(data,nboot,{bootfun,bootfun_args})
%  stats = bootknife(data,nboot,bootfun,alpha)
%  stats = bootknife(data,nboot,bootfun,alpha,strata)
%  stats = bootknife(data,nboot,bootfun,alpha,strata,nproc)
%  stats = bootknife(data,[2000,0],@mean,0.05,[],0)      % Default values
%  [stats,bootstat] = bootknife(...)
%  [stats,bootstat] = bootknife(...)
%  [stats,bootstat,bootsam] = bootknife(...)
%  bootknife(data,...);
%
%  stats = bootknife(data) resamples from the rows of a data sample (column 
%  vector or a matrix) and returns a structure with the following fields:
%    original: contains the result of applying bootfun to the data (x) 
%    bias: contains the bootstrap estimate of bias [7-8]
%    std_error: contains the bootstrap standard error
%    CI_lower: contains the lower bound of the bootstrap confidence interval
%    CI_upper: contains the upper bound of the bootstrap confidence interval
%  By default, the statistics relate to bootfun being @mean and the confidence
%  intervals are 95% bias-corrected and accelerated (BCa) intervals [1,4-5,9]. 
%
%  stats = bootknife(data,nboot) also specifies the number of bootstrap 
%  samples. nboot can be a scalar, or vector of upto two positive integers. 
%  By default, nboot is [2000,0], which implements a single bootstrap with 
%  the 2000 resamples, but larger numbers of resamples are recommended to  
%  reduce the Monte Carlo error, particularly for confidence intervals. If  
%  the second element of nboot is > 0, then the first and second elements  
%  of nboot correspond to the number of outer (first) and inner (second) 
%  bootstrap resamples respectively. Double bootstrap is used to improve 
%  the accuracy of the bias and the confidence intervals. For confidence 
%  intervals, this is achieved by calibrating the lower and upper interval 
%  ends to have tail probabilities of 2.5% and 97.5% [5]. Note that one 
%  can get away with a lower number of resamples in the second bootstrap 
%  to reduce the computational expense of the double bootstrap (e.g. [2000,
%  200]), since the algorithm uses linear interpolation to achieve near-
%  asymptotic calibration of confidence intervals [3]. The confidence 
%  intervals calculated (with either single or double bootstrap) are 
%  transformation invariant and have more accuracy and correctness 
%  compared to intervals derived from normal theory or to simple percentile
%  bootstrap confidence intervals.
%
%  stats = bootknife(data,nboot,bootfun) also specifies bootfun, a function 
%  handle, a string indicating the name of the function to apply to the data
%  (and each bootstrap resample), or a cell array where the first cell is the 
%  function handle or string, and other cells being arguments for that function, 
%  where the function must take data for the first input argument. bootfun can 
%  return a scalar value or vector. The default value(s) of bootfun is/are the
%  (column) mean(s). When bootfun is @mean or 'mean', residual narrowness
%  bias of central coverage is almost eliminated by using the Student's 
%  t-distribution to expand the percentiles before applying the BCa 
%  adjustments as described in [9].
%    Note that bootfun MUST calculate a statistic representative of the 
%  finite data sample, it should NOT be an estimate of a population 
%  parameter. For example, for the variance, set bootfun to {@var,1}, not 
%  @var or {@var,0}. Smooth functions of the data are preferable, (e.g. use
%  smoothmedian function instead of ordinary median). 
%    If single bootstrap is requested and bootfun cannot be executed during
%  leave-one-out jackknife, the acceleration constant will be set to 0 and
%  intervals will only be bias corrected. This warning can be disabled with the
%  following function call: warning('off','bootknife:jackfail'). 
%
%  stats = bootknife(data,nboot,bootfun,alpha) where alpha sets the lower 
%  and upper confidence interval ends. The value(s) in alpha must be between 0
%  and 1. If alpha is a scalar value, the nominal lower and upper percentiles
%  of the confidence are 100*(alpha/2)% and 100*(1-alpha/2)% respectively, and
%  nominal central coverage of the intervals is 100*(1-alpha)%. If alpha is a
%  vector with two elements, alpha becomes the quantiles for the confidence
%  intervals, and the intervals become percentile bootstrap confidence
%  intervals. If alpha is empty, NaN is returned for the confidence interval
%  ends. The default value for alpha is 0.05. 
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
%  [1] Hesterberg T.C. (2004) Unbiasing the Bootstrap—Bootknife Sampling 
%        vs. Smoothing; Proceedings of the Section on Statistics & the 
%        Environment. Alexandria, VA: American Statistical Association.
%  [2] Davison et al. (1986) Efficient Bootstrap Simulation.
%        Biometrika, 73: 555-66
%  [3] Gleason, J.R. (1988) Algorithms for Balanced Bootstrap Simulations. 
%        The American Statistician. Vol. 42, No. 4 pp. 263-266
%  [4] Efron (1987) Better Bootstrap Confidence Intervals. JASA, 
%        82(397): 171-185 
%  [5] Efron, and Tibshirani (1993) An Introduction to the
%        Bootstrap. New York, NY: Chapman & Hall
%  [6] Hall, Lee and Young (2000) Importance of interpolation when
%        constructing double-bootstrap confidence intervals. Journal
%        of the Royal Statistical Society. Series B. 62(3): 479-491
%  [7] Ouysee, R. (2011) Computationally efficient approximation for 
%        the double bootstrap mean bias correction. Economics Bulletin, 
%        AccessEcon, vol. 31(3), pages 2388-2403.
%  [8] Davison A.C. and Hinkley D.V (1997) Bootstrap Methods And Their 
%        Application. Chapter 3, pg. 104
%  [9] Hesterberg, Tim (2014), What Teachers Should Know about the 
%        Bootstrap: Resampling in the Undergraduate Statistics Curriculum, 
%        http://arxiv.org/abs/1411.5279
%
%  bootknife v1.9.0.0 (27/09/2022)
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
    nboot = [2000, 0];
  else
    if isempty(nboot)
      nboot = [2000, 0];
    else
      if ~isa (nboot, 'numeric')
        error('nboot must be numeric');
      end
      if any (nboot ~= abs (fix (nboot)))
        error ('nboot must contain positive integers')
      end    
      if (numel (nboot) > 2)
        error ('nboot cannot contain more than 2 values')
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
    if ~isa (alpha,'numeric') || numel (alpha) > 2
      error('alpha must be a numeric scalar value');
    end
    if any ((alpha < 0) | (alpha > 1))
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

  % Set number of outer and inner bootknife resamples
  B = nboot(1);
  if (numel (nboot) > 1)
    C =  nboot(2);
  else
    C = 0;
  end

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
    bootstat = zeros (m, B);
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
      l = []; % we do not have this information and it won't be the same for each column of data
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
        BOOTSAM = zeros (n, B, 'int32'); 
        for k = 1:K
          BOOTSAM(g(:, k),:) = boot (find (g(:, k)), B, true);
        end
      else
        % For more efficiency, if we don't need BOOTSAM, we can directly resample values of x
        BOOTSAM = [];
        X = zeros (n, B);
        for k = 1:K
          X(g(:, k),:) = boot (x(g(:, k),:), B, true);
        end
      end
    else
      if (nvar > 1) || (nargout > 2)
        % If we need BOOTSAM, can save some memory by making BOOTSAM an int32 datatype
        BOOTSAM = zeros (n, B, 'int32');
        BOOTSAM(:,:) = boot (n, B, true);
      else
        % For more efficiency, if we don't need BOOTSAM, we can directly resample values of x
        BOOTSAM = [];
        X = boot (x, B, true);
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
          bootstat = zeros (1, B);
          parfor b = 1:B; bootstat(b) = cellfunc (X(:, b)); end;
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
          bootstat = zeros (1, B);
          parfor b = 1:B; bootstat(b) = cellfunc (BOOTSAM(:, b)); end;
        end
      else
        % Evaluate bootfun on each bootstrap resample in SERIAL
        cellfunc = @(BOOTSAM) bootfun (x(BOOTSAM, :));
        bootstat = cellfun (cellfunc, num2cell (BOOTSAM, 1));
      end
    end
  end

  % Calculate the bootstrap bias, standard error and confidence intervals 
  if C > 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%% DOUBLE BOOTSTRAP %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (ncpus > 1)
      % PARALLEL execution of inner layer resampling for double (i.e. iterated) bootstrap
      if ISOCTAVE
        % OCTAVE
        % Set unique random seed for each parallel thread
        pararrayfun(ncpus, @boot, 1, 1, false, 1, 1:ncpus);
        if vectorized && isempty(BOOTSAM)
          cellfunc = @(x) bootknife (x, C, bootfun, [], strata, 0, T0, ISOCTAVE);
          bootout = parcellfun (ncpus, cellfunc, num2cell (X,1));
        else
          cellfunc = @(BOOTSAM) bootknife (x(BOOTSAM,:), C, bootfun, [], strata, 0, T0, ISOCTAVE);
          bootout = parcellfun (ncpus, cellfunc, num2cell (BOOTSAM,1));
        end
      else
        % MATLAB
        % Set unique random seed for each parallel thread
        parfor i = 1:ncpus; boot (1, 1, false, 1, i); end;
        % Perform inner layer of resampling
        bootout = struct ('original',zeros(1,B),...
                          'bias',zeros(1,B),...
                          'std_error',zeros(1,B),...
                          'CI_lower',zeros(1,B),...
                          'CI_upper',zeros(1,B),...
                          'Pr',zeros(1,B));
        if vectorized && isempty(BOOTSAM)
          cellfunc = @(x) bootknife (x, C, bootfun, [], strata, 0, T0, ISOCTAVE);
          parfor b = 1:B; bootout(b) = cellfunc (X(:,b)); end;
        else
          cellfunc = @(BOOTSAM) bootknife (x(BOOTSAM,:), C, bootfun, [], strata, 0, T0, ISOCTAVE);
          parfor b = 1:B; bootout(b) = cellfunc (BOOTSAM(:,b)); end;
        end
      end
    else
      % SERIAL execution of inner layer resampling for double bootstrap
      if vectorized && isempty(BOOTSAM)
        cellfunc = @(x) bootknife (x, C, bootfun, [], strata, 0, T0, ISOCTAVE);
        bootout = cellfun (cellfunc, num2cell (X,1));
      else
        cellfunc = @(BOOTSAM) bootknife (x(BOOTSAM,:), C, bootfun, [], strata, 0, T0, ISOCTAVE);
        bootout = cellfun (cellfunc, num2cell (BOOTSAM,1));
      end
    end
    Pr = cell2mat(arrayfun(@(S) S.Pr, bootout, 'UniformOutput', false));
    mu = cell2mat(arrayfun(@(S) S.bias, bootout, 'UniformOutput', false)) + ...
         cell2mat(arrayfun(@(S) S.original, bootout, 'UniformOutput', false));
    % Double bootstrap bias estimation
    b = mean (bootstat) - T0;
    c = mean (mu) - 2 * mean (bootstat) + T0;
    bias = b - c;
    % Bootstrap standard error
    se = std (bootstat, 1);
    if ~isempty(alpha)
      % Calibrate tail probabilities
      [cdf, u] = empcdf (Pr, 1);
      switch (numel (alpha))
        case 1
          % alpha is a two-tailed probability
          l = arrayfun ( @(p) interp1 (cdf, u, p, 'linear'), [alpha / 2, 1 - alpha / 2]);
        case 2
          % alpha is a vector of quantiles
          % Make sure quantiles are in the correct order
          if alpha(1) > alpha(2) 
            alpha = fliplr(alpha(:)');
          end
          l = arrayfun ( @(p) interp1 (cdf, u, p, 'linear'), alpha);
      end
      % Calibrated percentile bootstrap confidence intervals
      [cdf, t1] = empcdf (bootstat, 1);
      ci = arrayfun ( @(p) interp1 (cdf, t1, p, 'linear'), l);
    else
      ci = nan (1, 2);
    end
  else
    %%%%%%%%%%%%%%%%%%%%%%%%%%% SINGLE BOOTSTRAP %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bootstrap bias estimation
    bias = mean (bootstat) - T0;
    % Bootstrap standard error
    se = std (bootstat, 1);
    if ~isempty(alpha)
      switch (numel (alpha))
        case 1
          % Bias-corrected and accelerated bootstrap confidence intervals 
          % Create distribution functions
          stdnormcdf = @(x) 0.5 * (1 + erf (x / sqrt (2)));
          stdnorminv = @(p) sqrt (2) * erfinv (2 * p-1);
          % Calculate the median bias correction z0
          z0 = stdnorminv (sum (bootstat < T0)/ B);
          if isinf (z0) || isnan (z0)
            error('unable to calculate the bias correction z0')
          end
          % Use the Jackknife to calculate the acceleration constant
          try
            jackfun = @(i) bootfun (x(1:n ~= i, :));
            if (ncpus > 1)  
              % PARALLEL evaluation of bootfun on each jackknife resample 
              if ISOCTAVE
                % OCTAVE
                T = pararrayfun (ncpus, jackfun, 1:n);
              else
                % MATLAB
                T = zeros (n, 1);
                parfor i = 1:n; T(i) = jackfun (i); end;
              end
            else
              % SERIAL evaluation of bootfun on each jackknife resample
              T = arrayfun (jackfun, 1:n);
            end
            % Calculate empirical influence function
            if ~isempty(strata)
              gk = sum (g .* repmat (sum (g), n, 1), 2).';
              U = (gk - 1) .* (mean (T) - T);   
            else
              U = (n - 1) * (mean (T) - T);     
            end
            a = sum (U.^3) / (6 * sum (U.^2) ^ 1.5);
          catch
            a = 0;
            warning ('bootknife:jackfail','bootfun failed during jackknife, acceleration constant set to 0\n');
          end
          if strcmp (func2str (bootfun), 'mean')
            % If bootfun is the mean, expand percentiles using Student's 
            % t-distribution to improve central coverage for small samples
            if exist('betaincinv','file')
              studinv = @(p, df) - sqrt ( df ./ betaincinv (2 * p, df / 2, 0.5) - df);
            else
              % Earlier versions of matlab do not have betaincinv
              studinv = @(p, df) - sqrt ( df ./ betainv (2 * p, df / 2, 0.5) - df);
            end
            adj_alpha = stdnormcdf (studinv (alpha / 2, n - 1)) * 2;   
          else
            adj_alpha = alpha;
          end
          % Calculate BCa percentiles
          z1 = stdnorminv(adj_alpha / 2);
          z2 = stdnorminv(1 - adj_alpha / 2);
          l = cat (2, stdnormcdf (z0 + ((z0 + z1) / (1 - a * (z0 + z1)))),... 
                      stdnormcdf (z0 + ((z0 + z2) / (1 - a * (z0 + z2)))));
          [cdf, t1] = empcdf (bootstat, 1);
          ci = arrayfun ( @(p) interp1 (cdf, t1, p, 'linear'), l);
        case 2
          % alpha is a vector of quantiles
          % Make sure quantiles are in the correct order
          if alpha(1) > alpha(2) 
            alpha = fliplr(alpha(:)');
          end
          l = alpha;
          % Percentile bootstrap confidence intervals 
          [cdf, t1] = empcdf (bootstat, 1);
          ci = arrayfun ( @(p) interp1 (cdf, t1, p, 'linear'), l);
      end
    else
      ci = nan (1, 2);
    end
  end
  ci(alpha == 0) = -inf;
  ci(alpha == 1) = +inf;

  
  % Use quick interpolation to find the proportion (Pr) of bootstat <= REF
  if (nargin < 7)
    Pr = NaN;
  else
    if isempty(REF)
      Pr = NaN;
    else
      I = (bootstat <= REF);
      pr = sum (I);
      t = [max([min(bootstat), max(bootstat(I))]),...
           min([max(bootstat), min(bootstat(~I))])];
      if (pr < B) && ((t(2) - t(1)) > 0)
        % Linear interpolation to calculate Pr, which is required to calibrate alpha and improving confidence interval coverage 
        Pr = ((t(2) - REF) * pr / B + (REF - t(1)) * min ((pr + 1) / B, 1)) / (t(2) - t(1));
      else
        Pr = pr / B;
      end
    end
  end

  % Prepare stats output argument
  stats = struct;
  stats.original = T0;
  stats.bias = bias;
  stats.std_error = se;
  stats.CI_lower = ci(1);
  stats.CI_upper = ci(2);
  if ~isnan(Pr)
    stats.Pr = Pr;
  end
  
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
      fprintf (' Number of resamples (outer): %u \n', B);
      fprintf (' Number of resamples (inner): %u \n', C);
      if ~isempty(alpha)
        if (C > 0)
          fprintf (' Confidence interval type: Calibrated \n');
        else
          if (numel (alpha) > 1)
            fprintf (' Confidence interval type: Percentile \n');
          else
            fprintf (' Confidence interval type: Bias-corrected and accelerated (BCa) \n');
          end
        end
        if (numel (alpha) > 1)
          % alpha is a vector of quantiles
          coverage = 100*abs(alpha(2)-alpha(1));
        else
          % alpha is a two-tailed probability
          coverage = 100*(1-alpha);
        end
        if isempty (l)
          fprintf (' Confidence interval coverage: %g%%\n\n',coverage);
        else
          fprintf (' Confidence interval coverage: %g%% (%.1f%%, %.1f%%)\n\n',coverage,100*l);
        end
      end
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


%!test
%! ## Spatial Test Data from Table 14.1 of Efron and Tibshirani (1993)
%! ## An Introduction to the Bootstrap in Monographs on Statistics and Applied 
%! ## Probability 57 (Springer)
%! A = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!      0 33 28 34 4 32 24 47 41 24 26 30 41]';
%! ## Nonparametric 90% percentile confidence intervals (single bootstrap)
%! ## Table 14.2 percentile intervals are 100.8 - 233.9
%! boot (1, 1, true, [], 1); # Set random seed
%! stats = bootknife(A,2000,{@var,1},[0.05 0.95]);
%! assert (stats.original, 171.534023668639, 1e-09);
%! assert (stats.bias, -6.978681952662356, 1e-09);
%! assert (stats.std_error, 43.09155390135619, 1e-09);
%! assert (stats.CI_lower, 95.19578402366864, 1e-09);
%! assert (stats.CI_upper, 238.9609467455621, 1e-09);
%! ## Nonparametric 90% BCa confidence intervals (single bootstrap)
%! ## Table 14.2 BCa intervals are 115.8 - 259.6
%! boot (1, 1, true, [], 1); # Set random seed
%! stats = bootknife(A,2000,{@var,1},0.1);
%! assert (stats.original, 171.534023668639, 1e-09);
%! assert (stats.bias, -6.978681952662356, 1e-09);
%! assert (stats.std_error, 43.09155390135619, 1e-09);
%! assert (stats.CI_lower, 115.6455796312253, 1e-09);
%! assert (stats.CI_upper, 269.4469269661803, 1e-09);
%! ## Nonparametric 90% calibrated confidence intervals (double bootstrap)
%! boot (1, 1, true, [], 1); # Set random seed
%! stats = bootknife(A,[2000,200],{@var,1},0.1);
%! assert (stats.original, 171.534023668639, 1e-09);
%! assert (stats.bias, -7.407638883135036, 1e-09);
%! assert (stats.std_error, 43.09155390135619, 1e-09);
%! assert (stats.CI_lower, 111.39427003007, 1e-09);
%! assert (stats.CI_upper, 313.7115384615385, 1e-09);
%! # Exact intervals based on theory are 118.4 - 305.2 (Table 14.2)
%! # Note that all of the nootknife intervals are slightly wider than the
%! # non-parametric intervals in Table 14.2 because the bootknife (rather than
%! # standard bootstrap) resampling used here reduces small sample bias

%!test
%! ## Data from Table 14.1: Spatial Test Data in DiCiccio and Efron (1996)
%! ## Bootstrap Confidence Intervals. Statistical Science. 11(3):189-228
%! baseline = [2.12,4.35,3.39,2.51,4.04,5.1,3.77,3.35,4.1,3.35, ...
%!             4.15,3.56, 3.39,1.88,2.56,2.96,2.49,3.03,2.66,3]';
%! oneyear  = [2.47,4.61,5.26,3.02,6.36,5.93,3.93,4.09,4.88,3.81, ...
%!             4.74,3.29,5.55,2.82,4.23,3.23,2.56,4.31,4.37,2.4]';
%! bootfun = @(M) corr(M(:,1),M(:,2));
%! ## Nonparametric 90% BCa confidence intervals (single bootstrap)
%! ## Table 2 BCa intervals are 0.55 - 0.85
%! boot (1, 1, true, [], 1); # Set random seed
%! stats = bootknife([baseline,oneyear],[2000,0],bootfun,0.1);
%! assert (stats.original, 0.7231653678920302, 1e-09);
%! assert (stats.bias, -0.007860430929598206, 1e-09);
%! assert (stats.std_error, 0.09328837054352047, 1e-09);
%! assert (stats.CI_lower, 0.5477225147834641, 1e-09);
%! assert (stats.CI_upper, 0.8457573378934136, 1e-09);
%! # Exact intervals based on theory are 0.47 - 0.86 (Table 14.2)