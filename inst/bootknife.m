%  Function File: bootknife
%
%  Bootknife (bootstrap) resampling and statistics
%
%  This function takes a DATA sample (containing n rows) and uses bootstrap
%  methods to calculate a bias of the parameter estimate, a standard error,
%  and 95% confidence intervals. Specifically, the method uses bootknife
%  resampling [1], which involves creating leave-one-out jackknife samples
%  of size n - 1 and then drawing samples of size n with replacement from the
%  jackknife samples. The resampling of DATA rows is balanced in order to
%  reduce Monte Carlo error [2,3] and the confidence intervals are obtained
%  from the quantiles of a kernel density estimate of bootstrap statistics (with
%  shrinkage corrrection). By default, the confidence intervals are bias-
%  corrected and accelerated (BCa) [4-5]. BCa intervals are fast to compute and
%  have good coverage and correctness when combined with bootknife resampling
%  as it is here [1]. If double bootstrap is requested, the algorithm uses
%  calibration to improve the accuracy of the bias and standard error, and the
%  coverage of the confidence intervals [6-12]. 
%
%  STATS = bootknife (DATA)
%  STATS = bootknife (DATA, NBOOT)
%  STATS = bootknife (DATA, NBOOT, BOOTFUN)
%  STATS = bootknife (DATA, NBOOT, {BOOTFUN, ...})
%  STATS = bootknife (DATA, NBOOT, ..., ALPHA)
%  STATS = bootknife (DATA, NBOOT, ..., ALPHA, STRATA)
%  STATS = bootknife (DATA, NBOOT, ..., ALPHA, STRATA, NPROC)
%  [STATS, BOOTSTAT] = bootknife (...)
%  [STATS, BOOTSTAT] = bootknife (...)
%  [STATS, BOOTSTAT, BOOTSAM] = bootknife (...)
%  bootknife (DATA,...);
%  bootknife (DATA, [2000, 0], @mean, 0.05, [], 0)            % Defaults (single)
%  bootknife (DATA, [2000, 200], @mean, [0.025,0.975], [], 0) % Defaults (double)
%
%  STATS = bootknife (DATA) resamples from the rows of a DATA sample (column 
%  vector or a matrix) and returns a structure with the following fields:
%    original: contains the result of applying BOOTFUN to the DATA 
%    bias: contains the bootstrap estimate of bias [11-12]
%    std_error: contains the bootstrap standard error
%    CI_lower: contains the lower bound of the bootstrap confidence interval
%    CI_upper: contains the upper bound of the bootstrap confidence interval
%  By default, the statistics relate to BOOTFUN being @mean and the confidence
%  intervals are 95% bias-corrected and accelerated (BCa) intervals [1,4-5].
%  If DATA is a cell array of column vectors, the vectors are passed to BOOTFUN
%  as separate input arguments.
%
%  STATS = bootknife (DATA, NBOOT) also specifies the number of bootstrap 
%  samples. NBOOT can be a scalar, or vector of upto two positive integers. 
%  By default, NBOOT is [2000,0], which implements a single bootstrap with 
%  the 2000 resamples, but larger numbers of resamples are recommended to  
%  reduce the Monte Carlo error, particularly for confidence intervals. If  
%  the second element of NBOOT is > 0, then the first and second elements  
%  of NBOOT correspond to the number of outer (first) and inner (second) 
%  bootstrap resamples respectively. This so called double bootstrap is used
%  the accuracy of the bias, standard error and confidence intervals. The
%  latter is achieved by calibrating the lower and upper interval ends to
%  have nominal tail probabilities of 2.5% and 97.5% [5]. Note that one 
%  can get away with a lower number of resamples in the second bootstrap 
%  to reduce the computational expense of the double bootstrap (e.g. [2000,
%  200]), since the algorithm uses linear interpolation to achieve near-
%  asymptotic calibration of confidence intervals [3]. The confidence 
%  intervals calculated (with either single or double bootstrap) are 
%  transformation invariant and can have better coverage and/or accuracy
%  compared to intervals derived from normal theory or to simple percentile
%  bootstrap confidence intervals [5-10].
%
%  STATS = bootknife (DATA, NBOOT, BOOTFUN) also specifies BOOTFUN, a function 
%  handle, a string indicating the name of the function to apply to the DATA
%  (and each bootstrap resample), or a cell array where the first cell is the 
%  function handle or string, and other cells being additional input arguments 
%  for BOOTFUN, where BOOTFUN must take DATA for the first input arguments.
%  BOOTFUN can return a scalar value or vector. BOOTFUN must calculate a
%  statistic representative of the finite DATA sample, it should not be an
%  unbiased estimate of a population parameter. For example, for the variance,
%  set BOOTFUN to {@var,1} (not @var or {@var,0}), and for the correlation
%  coefficient, use the provided @cor function. Smooth functions of the DATA
%  are preferable, (e.g. use smoothmedian function instead of the ordinary
%  median). The default value(s) of BOOTFUN is/are the (column) mean(s).
%    When BOOTFUN is @mean or 'mean', residual narrowness bias of central
%  coverage is almost eliminated by using Student's t-distribution to expand
%  the nominal tail probabilities [13].
%    Note that if single bootstrap is requested and BOOTFUN cannot be executed
%  during leave-one-out jackknife, the acceleration constant will be set to 0
%  and intervals will be bias-corrected only.
%
%  STATS = bootknife (DATA, NBOOT, BOOTFUN, ALPHA) where ALPHA is numeric and
%  sets the lower and upper bounds of the confidence interval(s). The value(s)
%  of ALPHA must be between 0 and 1. ALPHA can either be:
%
%  1) a scalar value to set the (nominal) central coverage (e.g. .05) to
%  100*(1-ALPHA)% with (nominal) lower and upper percentiles of the confidence
%  intervals at 100*(ALPHA/2)% and 100*(1-ALPHA/2)% respectively, or
%
%  2) a vector containing a pair of probabilities (e.g. [.025, .975]) to set
%  the (nominal) lower and upper percentiles of the confidence interval(s) at
%  100*(ALPHA(1))% and 100*(ALPHA(2))%.
%
%  The method for constructing confidence intervals is determined by the
%  combined settings of ALPHA and NBOOT:
%
%  - PERCENTILE: ALPHA must be a pair of probabilities and NBOOT must be a
%    scalar value (or the second element in NBOOT must be zero).
%
%  - BIAS-CORRECTED AND ACCELERATED (BCa): ALPHA must be a scalar value and
%    NBOOT must be a scalar value (or the second element in NBOOT must be zero).
%
%  - CALIBRATED PERCENTILE (symmetric, equal-tailed): ALPHA must be a scalar
%    value and NBOOT must be a vector of two positive, non-zero integers (for
%    double bootstrap). The method used corresponds to the 2-sided intervals in
%    [7,9].
%
%  - CALIBRATED PERCENTILE (asymmetric): ALPHA must be must be a pair of
%    probabilities and NBOOT must be a vector of two positive, non-zero integers
%    (for double bootstrap). The method used corresponds to the 1-sided (lower
%    and upper) intervals in [7-10], which are a fast and accurate approximation
%    of the confidence point calibration algorithm 18.1 in [5].
%
%  Confidence interval endpoints are not calculated when the value(s) of ALPHA
%  is/are NaN. If empty (or not specified), the default value of ALPHA is 0.05
%  for single bootstrap (i.e. 95% BCa intervals), or [0.025, 0.975] for double 
%  bootstrap (i.e. 95% calibrated intervals, with 2.5% and 97.5% confidence
%  point calibration). 
%
%  STATS = bootknife (DATA, NBOOT, BOOTFUN, ALPHA, STRATA) also sets STRATA, 
%  which are identifiers that define the grouping of the DATA rows
%  for stratified bootstrap resampling. STRATA should be a column vector 
%  or cell array the same number of rows as the DATA. When resampling is 
%  stratified, the groups (or stata) of DATA are equally represented across 
%  the bootstrap resamples. If this input argument is not specified or is 
%  empty, no stratification of resampling is performed. 
%
%  STATS = bootknife (DATA, NBOOT, BOOTFUN, ALPHA, STRATA, NPROC) sets the
%  number of parallel processes to use to accelerate computations on 
%  multicore machines, specifically non-vectorized function evaluations,
%  double bootstrap resampling and jackknife function evaluations. This
%  feature requires the Parallel package (in Octave), or the Parallel
%  Computing Toolbox (in Matlab).
%
%  [STATS, BOOTSTAT] = bootknife (...) also returns BOOTSTAT, a vector of
%  bootstrap statistics calculated over the (first, or outer layer of)
%  bootstrap resamples.
%
%  [STATS, BOOTSTAT, BOOTSAM] = bootknife (...) also returns BOOTSAM, the
%  matrix of indices (32-bit integers) used for the (first, or outer
%  layer of) bootstrap resampling. Each column in BOOTSAM corresponds
%  to one bootstrap resample and contains the row indices of the values
%  drawn from the nonscalar DATA argument to create that sample.
%
%  bootknife (DATA, ...); returns a pretty table of the output including
%  the bootstrap settings and the result of evaluating BOOTFUN on the
%  DATA along with bootstrap estimates of bias, standard error, and
%  lower and upper 100*(1-ALPHA)% confidence limits.
%
%  Requirements: The function file boot.m (or better boot.mex) also
%  distributed in the statistics-bootstrap package.
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
%  [6] Beran (1987). Prepivoting to Reduce Level Error of Confidence Sets.
%        Biometrika, 74(3), 457–468.
%  [7] Lee and Young (􏰒1999) The e􏰚ffect of Monte Carlo approximation on coverage
%        error of double-bootstrap con®dence intervals. J R Statist Soc B. 
%        61:353-366.
%  [8] Booth J. and Presnell B. (1998) Allocation of Monte Carlo Resources for
%        the Iterated Bootstrap. J. Comput. Graph. Stat. 7(1):92-112 
%  [9] Hall, Lee and Young (2000) Importance of interpolation when
%        constructing double-bootstrap confidence intervals. Journal
%        of the Royal Statistical Society. Series B. 62(3): 479-491
%  [10] DiCiccio, Martin and Young (1992) Fast and accurate approximate double
%        bootstrap confidence intervals. Biometrika. 79(2):285-95
%  [11] Ouysee, R. (2011) Computationally efficient approximation for 
%        the double bootstrap mean bias correction. Economics Bulletin, 
%        AccessEcon, vol. 31(3), pages 2388-2403.
%  [12] Davison A.C. and Hinkley D.V (1997) Bootstrap Methods And Their 
%        Application. Chapter 3, pg. 104
%  [13] Hesterberg, Tim (2014), What Teachers Should Know about the 
%        Bootstrap: Resampling in the Undergraduate Statistics Curriculum, 
%        http://arxiv.org/abs/1411.5279
%
%  bootknife (version 2022.12.25)
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


function [stats, bootstat, BOOTSAM] = bootknife (x, nboot, bootfun, alpha, ...
                              strata, ncpus, REF, ISOCTAVE, BOOTSAM, ERRCHK)

  % Input argument names in all-caps are for internal use only
  % REF, ISOCTAVE, BOOTSAM, ERRCHK are undocumented input arguments required
  % for some of the functionalities of bootknife

  % Store local functions in a stucture for parallel processes
  localfunc = struct ('col2args', @col2args, ...
                      'empcdf', @empcdf, ...
                      'kdeinv', @kdeinv, ...
                      'ExpandProbs', @ExpandProbs);

  % Set defaults and check for errors (if applicable)
  if ((nargin < 10) || isempty (ERRCHK) || ERRCHK)
    if (nargin < 1)
      error ('bootknife: DATA must be provided');
    end
    if ((nargin < 2) || isempty(nboot))
      nboot = [2000, 0];
    else
      if (~ isa (nboot, 'numeric'))
        error ('bootknife: NBOOT must be numeric');
      end
      if (numel (nboot) > 2)
        error ('bootknife: NBOOT cannot contain more than 2 values');
      end
      if (any (nboot ~= abs (fix (nboot))))
        error ('bootknife: NBOOT must contain positive integers');
      end    
      if (numel(nboot) == 1)
        nboot = [nboot, 0];
      end
    end
    if ((nargin < 3) || isempty (bootfun))
      bootfun = @mean;
    else
      if (iscell (bootfun))
        func = bootfun{1};
        args = bootfun(2:end);
        bootfun = @(varargin) func (varargin{:}, args{:});
      end
      if (ischar (bootfun))
        % Convert character string of a function name to a function handle
        bootfun = str2func (bootfun);
      end
      if (~ isa (bootfun, 'function_handle'))
        error ('bootknife: BOOTFUN must be a function name or function handle');
      end
    end
    % Store bootfun as string for printing output at the end
    bootfun_str = func2str (bootfun);
    if (iscell (x))
      % If DATA is a cell array of equal size colunmn vectors, convert the cell
      % array to a matrix and redefine bootfun to parse multiple input arguments
      x = [x{:}];
      nvar = size (x, 2);
      bootfun = @(x) localfunc.col2args (bootfun, x, nvar);
    end
    if (~ (size (x, 1) > 1))
      error ('bootknife: DATA must contain more than one row');
    end
    if ((nargin < 4) || isempty (alpha))
      if (nboot(2) > 0)
        alpha = [0.025, 0.975];
        nalpha = 2;
      else
        alpha = 0.05;
        nalpha = 1;
      end
    else
      nalpha = numel (alpha);
      if (~isa (alpha,'numeric') || (nalpha > 2))
        error ('bootknife: ALPHA must be a scalar (two-tailed probability) or a vector (pair of probabilities)');
      end
      if any ((alpha < 0) | (alpha > 1))
        error ('bootknife: value(s) in ALPHA must be between 0 and 1');
      end
      if (nalpha > 1)
        % alpha is a pair of probabilities
        % Make sure probabilities are in the correct order
        if (alpha(1) > alpha(2) )
          error ('bootknife: the pair of probabilities must be in ascending numeric order');
        end
      end
    end
    if ((nargin < 5) || isempty (strata))
      strata = [];
    else  
      if (size (strata, 1) ~= size (x, 1))
        error ('bootknife: STRATA should be a column vector or cell array with the same number of rows as the DATA');
      end
    end
    if ((nargin < 6) || isempty (ncpus)) 
      ncpus = 0;    % Ignore parallel processing features
    else
      if ~isa (ncpus, 'numeric')
        error('bootknife: NPROC must be numeric');
      end
      if any (ncpus ~= abs (fix (ncpus)))
        error ('bootknife: NPROC must be a positive integer');
      end    
      if (numel (ncpus) > 1)
        error ('bootknife: NPROC must be a scalar value');
      end
    end
    if ((nargin < 8) || isempty (ISOCTAVE))
      % Check if running in Octave (else assume Matlab)
      info = ver; 
      ISOCTAVE = any (ismember ({info.Name}, 'Octave'));
    end
    if ISOCTAVE
      ncpus = min (ncpus, nproc);
    else
      ncpus = min (ncpus, feature ('numcores'));
    end
  end

  % Determine properties of the DATA (x)
  [n, nvar] = size (x);
  if (n < 2)
    error ('bootknife: the data must be a data vector or matrix with > 1 row')
  end

  % Set number of outer and inner bootknife resamples
  B = nboot(1);
  if (numel (nboot) > 1)
    C = nboot(2);
  else
    C = 0;
  end

  % Evaluate bootfun on the DATA
  T0 = bootfun (x);
  if all (size (T0) > 1)
    error ('bootknife: BOOTFUN must return either a scalar or a vector');
  end

  % Check whether bootfun is vectorized
  if (nvar > 1)
    M = cell2mat (cellfun (@(i) repmat (x(:,i), 1, 2), ...
                  num2cell (1:nvar), 'UniformOutput', false));
  else
    M = repmat (x, 1, 2);
  end
  try
    chk = bootfun (M);
    if (all (size (chk) == [1, 2]) && all (chk == bootfun (x)))
      vectorized = true;
    else
      vectorized = false;
    end
  catch
    vectorized = false;
  end

  % Initialize probabilities
  l = [];

  % If applicable, check we have parallel computing capabilities
  if (ncpus > 1)
    if ISOCTAVE
      pat = '^parallel';
      software = pkg('list');
      names = cellfun (@(S) S.name, software, 'UniformOutput', false);
      status = cellfun (@(S) S.loaded, software, 'UniformOutput', false);
      index = find (~cellfun (@isempty, regexpi (names,pat)));
      if (~ isempty (index))
        if logical (status{index})
          PARALLEL = true;
        else
          PARALLEL = false;
        end
      else
        PARALLEL = false;
      end
    else
      info = ver; 
      if ismember ('Parallel Computing Toolbox', {info.Name})
        PARALLEL = true;
      else
        PARALLEL = false;
      end
    end
  end
  
  % If applicable, setup a parallel pool (required for MATLAB)
  if (~ ISOCTAVE)
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
        warning ('MATLAB Parallel Computing Toolbox is not installed or operational. Falling back to serial processing.');
        ncpus = 1;
      end
    end
  else
    if ((ncpus > 1) && ~ PARALLEL)
      if ISOCTAVE
        % OCTAVE Parallel Computing Package is not installed or loaded
        warning ('OCTAVE Parallel Computing Package is not installed and/or loaded. Falling back to serial processing.');
      else
        % MATLAB Parallel Computing Toolbox is not installed or loaded
        warning ('MATLAB Parallel Computing Toolbox is not installed and/or loaded. Falling back to serial processing.');
      end
      ncpus = 0;
    end
  end

  % If the function of the DATA is a matrix, calculate and return the bootstrap statistics for each column 
  sz = size (T0);
  m = prod (sz);
  if (m > 1)
    if (nvar == m)
      % Confirm whether bootfun is vectorized when bootfun is evaluated across
      % all columns of x, considering each of them as univariate DATA vectors
      chk = bootfun (cat (2, x(:,1), x(:,1)));
      if ( all (size (chk) == [1, 2]) && all (chk == bootfun (x(:,1))) )
        vectorized = true;
      else
        vectorized = false;
      end
    end
    % Use bootknife for each element of the output of bootfun
    % Note that row indices in the resamples are the same for all columns of DATA
    stats = struct ('original', zeros(sz),...
                    'bias', zeros(sz),...
                    'std_error', zeros(sz),...
                    'CI_lower', zeros(sz),...
                    'CI_upper', zeros(sz));
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
      if (nalpha > 1) && (C == 0)
        l = alpha;
      end
      print_output (stats, B, C, alpha, l, m, bootfun_str);
    else
      [warnmsg, warnID] = lastwarn;
      if ismember (warnID, {'bootknife:biasfail','bootknife:jackfail'})
        warning ('bootknife:lastwarn', warnmsg);
      end
      lastwarn ('', '');
    end
    return
  end

  % Evaluate strata input argument
  if (~ isempty (strata))
    if ~isnumeric (strata)
      % Convert strata to numeric ID
      [jnk1, jnk2, strata] = unique (strata);
      clear jnk1 jnk2;
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
    K = 1;
  end

  % Perform balanced bootknife resampling
  unbiased = true;  % Set to true for bootknife resampling
  if ((nargin < 9) || isempty (BOOTSAM))
    if ~isempty (strata)
      if (nvar > 1) || (nargout > 2)
        % If we need BOOTSAM, can save some memory by making BOOTSAM an int32 datatype
        BOOTSAM = zeros (n, B, 'int32'); 
        for k = 1:K
          BOOTSAM(g(:, k),:) = boot (find (g(:, k)), B, unbiased);
        end
      else
        % For more efficiency, if we don't need BOOTSAM, we can directly resample values of x
        BOOTSAM = [];
        X = zeros (n, B);
        for k = 1:K
          X(g(:, k),:) = boot (x(g(:, k),:), B, unbiased);
        end
      end
    else
      if (nvar > 1) || (nargout > 2)
        % If we need BOOTSAM, can save some memory by making BOOTSAM an int32 datatype
        BOOTSAM = zeros (n, B, 'int32');
        BOOTSAM(:,:) = boot (n, B, unbiased);
      else
        % For more efficiency, if we don't need BOOTSAM, we can directly resample values of x
        BOOTSAM = [];
        X = boot (x, B, unbiased);
      end
    end
  end

  % Evaluate bootfun each bootstrap resample
  if (isempty (BOOTSAM))
    if vectorized
      % Vectorized evaluation of bootfun on the DATA resamples
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
          parfor b = 1:B; bootstat(b) = bootfun (X(:, b)); end
        end
      else
        bootstat = cellfun (bootfun, num2cell (X, 1));
      end
    end
  else
    if vectorized
      % DATA resampling (using BOOTSAM) and vectorized evaluation of bootfun on 
      % the DATA resamples 
      if (nvar > 1)
        % Multivariate
        % Perform DATA sampling
        X = cell2mat (cellfun (@(i) reshape (x(BOOTSAM,i), n, B), ...
                      num2cell (1:nvar), 'UniformOutput', false));
      else
        % Univariate
        % Perform DATA sampling
        X = x(BOOTSAM);
      end
        % Function evaluation on bootknife samples
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
          parfor b = 1:B; bootstat(b) = cellfunc (BOOTSAM(:, b)); end
        end
      else
        % Evaluate bootfun on each bootstrap resample in SERIAL
        bootstat = cellfun (cellfunc, num2cell (BOOTSAM, 1));
      end
    end
  end

  % Remove bootstrap statistics that are NaN and their associated DATA resamples
  ridx = isnan (bootstat);
  bootstat(ridx) = [];
  if isempty(BOOTSAM)
    X(:,ridx) = [];
  else
    BOOTSAM(:,ridx) = [];
  end
  B = B - sum (ridx);

  % Calculate the bootstrap bias, standard error and confidence intervals 
  if (C > 0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%% DOUBLE BOOTSTRAP %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (ncpus > 1)
      % PARALLEL execution of inner layer resampling for double (i.e. iterated) bootstrap
      if ISOCTAVE
        % OCTAVE
        % Set unique random seed for each parallel thread
        pararrayfun (ncpus, @boot, 1, 1, false, 1:ncpus);
        if vectorized && isempty(BOOTSAM)
          cellfunc = @(x) bootknife (x, C, bootfun, NaN, strata, 0, T0, ISOCTAVE, [], false);
          bootout = parcellfun (ncpus, cellfunc, num2cell (X,1));
        else
          cellfunc = @(BOOTSAM) bootknife (x(BOOTSAM,:), C, bootfun, NaN, strata, 0, T0, ISOCTAVE, [], false);
          bootout = parcellfun (ncpus, cellfunc, num2cell (BOOTSAM,1));
        end
      else
        % MATLAB
        % Set unique random seed for each parallel thread
        parfor i = 1:ncpus; boot (1, 1, false, i); end
        % Perform inner layer of resampling
        bootout = struct ('original', zeros (1,B),...
                          'bias', zeros (1,B),...
                          'std_error', zeros (1,B),...
                          'CI_lower', zeros (1,B),...
                          'CI_upper', zeros (1,B),...
                          'Pr', zeros (1,B));
        if vectorized && isempty(BOOTSAM)
          cellfunc = @(x) bootknife (x, C, bootfun, NaN, strata, 0, T0, ISOCTAVE, [], false);
          parfor b = 1:B; bootout(b) = cellfunc (X(:,b)); end
        else
          cellfunc = @(BOOTSAM) bootknife (x(BOOTSAM,:), C, bootfun, NaN, strata, 0, T0, ISOCTAVE, [], false);
          parfor b = 1:B; bootout(b) = cellfunc (BOOTSAM(:,b)); end
        end
      end
    else
      % SERIAL execution of inner layer resampling for double bootstrap
      if vectorized && isempty(BOOTSAM)
        cellfunc = @(x) bootknife (x, C, bootfun, NaN, strata, 0, T0, ISOCTAVE, [], false);
        bootout = cellfun (cellfunc, num2cell (X,1));
      else
        cellfunc = @(BOOTSAM) bootknife (x(BOOTSAM,:), C, bootfun, NaN, strata, 0, T0, ISOCTAVE, [], false);
        bootout = cellfun (cellfunc, num2cell (BOOTSAM,1));
      end
    end
    % Double bootstrap bias estimation
    mu = cell2mat (arrayfun (@(S) S.bias, bootout, 'UniformOutput', false)) + ...
         cell2mat (arrayfun (@(S) S.original, bootout, 'UniformOutput', false));
    b = mean (bootstat) - T0;
    c = mean (mu) - 2 * mean (bootstat) + T0;
    bias = b - c;
    % Double bootstrap multiplicative correction of the standard error
    V = cell2mat (arrayfun (@(S) S.std_error^2, bootout, 'UniformOutput', false));
    se = sqrt (var (bootstat)^2 / mean (V));
    % Double bootstrap confidence intervals
    if (~ isnan (alpha))
      U = cell2mat (arrayfun (@(S) S.Pr, bootout, 'UniformOutput', false));
      % Calibrate tail probabilities
      switch (nalpha)
        case 1
          % alpha is a two-tailed probability (scalar)
          % Calibrate central coverage and construct equal-tailed intervals (2-sided)
          [cdf, v] = localfunc.empcdf (abs (2 * U - 1), 1);
          vk = interp1 (cdf, v, 1 - alpha, 'linear');
          l = arrayfun (@(sign) 0.5 * (1 + sign * vk), [-1, 1]);
        case 2
          % alpha is a pair of probabilities (vector)
          % Calibrate coverage but construct endpoints separately (1-sided)
          % This is equivalent to algorithm 18.1 in Efron, and Tibshirani (1993)
          [cdf, u] = localfunc.empcdf (U, 1);
          l = arrayfun ( @(p) interp1 (cdf, u, p, 'linear'), alpha);
      end
      % Calibrated percentile bootstrap confidence intervals
      % Kernel density estimate (with shrinkage correction)
      ci = localfunc.kdeinv (l, bootstat, se * sqrt (1 / (n - K)), 1 - 1 / (n - K));
      % Linear interpolation (legacy)
      %[cdf, t1] = localfunc.empcdf (bootstat, 1);
      %ci = arrayfun ( @(p) interp1 (cdf, t1, p, 'linear'), l);
    else
      ci = nan (1, 2);
    end
  else
    %%%%%%%%%%%%%%%%%%%%%%%%%%% SINGLE BOOTSTRAP %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bootstrap bias estimation
    bias = mean (bootstat) - T0;
    % Bootstrap standard error
    se = std (bootstat);  % Unbiased estimate since we used bootknife resampling
    if (~ isnan (alpha))
      % If bootfun is the arithmetic meam, perform expansion based on t-distribution
      if (strcmp (func2str (bootfun), 'mean'))
        alpha = (3 - nalpha) * localfunc.ExpandProbs (alpha / (3 - nalpha), n - K);
      end
      state = warning;
      if ISOCTAVE
        warning ('on', 'quiet');
      else
        warning off;
      end
      % Create distribution functions
      stdnormcdf = @(x) 0.5 * (1 + erf (x / sqrt (2)));
      stdnorminv = @(p) sqrt (2) * erfinv (2 * p-1);
      switch (nalpha)
        case 1
          % Attempt to form bias-corrected and accelerated (BCa) bootstrap confidence intervals. 
          % Use the Jackknife to calculate the acceleration constant (a)
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
                parfor i = 1:n; T(i) = jackfun (i); end
              end
            else
              % SERIAL evaluation of bootfun on each jackknife resample
              T = arrayfun (jackfun, 1:n);
            end
            % Calculate empirical influence function
            if (~ isempty (strata))
              gk = sum (g .* repmat (sum (g), n, 1), 2).';
              U = (gk - 1) .* (mean (T) - T);   
            else
              U = (n - 1) * (mean (T) - T);
            end
            a = sum (U.^3) / (6 * sum (U.^2) ^ 1.5);
          catch
            % Revert to bias-corrected (BC) bootstrap confidence intervals
            warning ('bootknife:jackfail','bootfun failed during jackknife, acceleration constant set to 0\n');
            a = 0;
          end
          % Calculate the bias correction constant (z0)
          % Calculate the median bias correction z0
          z0 = stdnorminv (sum (bootstat < T0) / B);
          if (isinf (z0) || isnan (z0))
            % Revert to percentile bootstrap confidence intervals
            warning ('bootknife:biasfail','unable to calculate the bias correction constant, reverting to percentile intervals\n');
            z0 = 0;
            a = 0; 
            l = [alpha / 2, 1 - alpha / 2];
          end
          if (isempty (l))
            % Calculate BCa or BC percentiles
            z1 = stdnorminv (alpha / 2);
            z2 = stdnorminv (1 - alpha / 2);
            l = cat (2, stdnormcdf (z0 + ((z0 + z1) / (1 - a * (z0 + z1)))),... 
                        stdnormcdf (z0 + ((z0 + z2) / (1 - a * (z0 + z2)))));
          end
        case 2
          % alpha is a vector of probabilities
          l = alpha;
      end
      % Intervals constructed from kernel density estimate of the bootstrap
      % (with shrinkage correction)
      ci = localfunc.kdeinv (l, bootstat, se * sqrt (1 / (n - K)), 1 - 1 / (n - K));
      warning (state);
      if ISOCTAVE
        warning ('off', 'quiet');
      end
    else
      ci = nan (1, 2);
    end
  end
  ci(l == 0) = -inf; 
  ci(l == 1) = +inf;

  % Use quick interpolation to find the proportion (Pr) of bootstat <= REF
  if ((nargin < 7) || isempty (REF))
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
        % Linear interpolation to calculate Pr, which is required to calibrate 
        % alpha and improving confidence interval coverage 
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
  if (~ isnan (Pr))
    stats.Pr = Pr;
  end
  
  % Print output if no output arguments are requested
  if (nargout == 0) 
    print_output (stats, B, C, alpha, l, m, bootfun_str);
  else
    if isempty(BOOTSAM)
      [warnmsg, warnID] = lastwarn;
      if ismember (warnID, {'bootknife:biasfail','bootknife:jackfail'})
        warning ('bootknife:lastwarn', warnmsg);
      end
      lastwarn ('', '');
    end
  end

end

%--------------------------------------------------------------------------

function print_output (stats, B, C, alpha, l, m, bootfun_str)

    fprintf (['\nSummary of non-parametric bootstrap estimates of bias and precision\n',...
              '******************************************************************************\n\n']);
    fprintf ('Bootstrap settings: \n');
    fprintf (' Function: %s\n', bootfun_str);
    fprintf (' Resampling method: Balanced, bootknife resampling \n');
    fprintf (' Number of resamples (outer): %u \n', B);
    fprintf (' Number of resamples (inner): %u \n', C);
    if (~ isempty (alpha) && ~ all (isnan (alpha)))
      nalpha = numel (alpha);
      if (C > 0)
        if (nalpha > 1)
          fprintf (' Confidence interval type: Calibrated percentile (asymmetric, 1-sided)\n');
        else
          fprintf (' Confidence interval type: Calibrated percentile (symmetric, 2-sided)\n');
        end
      else
        if (nalpha > 1)
          fprintf (' Confidence interval type: Percentile \n');
        else
          [jnk, warnID] = lastwarn;
          switch warnID
            case 'bootknife:biasfail'
              fprintf (' Confidence interval type: Percentile \n');
            case 'bootknife:jackfail'
              fprintf (' Confidence interval type: Bias-corrected (BC)\n');
            otherwise
              fprintf (' Confidence interval type: Bias-corrected and accelerated (BCa) \n');
          end
        end
      end
      if (nalpha > 1)
        % alpha is a vector of probabilities
        coverage = 100 * abs (alpha(2) - alpha(1));
      else
        % alpha is a two-tailed probability
        coverage = 100 * (1 - alpha);
      end
      if (isempty (l))
        fprintf (' Confidence interval coverage: %.3g%%\n\n',coverage);
      else
        fprintf (' Confidence interval coverage: %.3g%% (%.1f%%, %.1f%%)\n\n',coverage,100*l);
      end
    end
    fprintf ('Bootstrap Statistics: \n');
    fprintf (' original       bias           std_error      CI_lower       CI_upper    \n');
    for i = 1:m
      fprintf (' %#-+12.6g   %#-+12.6g   %#-+12.6g   %#-+12.6g   %#-+12.6g \n',... 
               [stats(i).original, stats(i).bias, stats(i).std_error, stats(i).CI_lower, stats(i).CI_upper]);
    end
    fprintf ('\n');
    lastwarn ('', '');  % reset last warning

end

%--------------------------------------------------------------------------

function retval = col2args (func, x, nvar)

  % Usage: retval = col2args (func, x, nvar)
  % col2args evaluates func on the columns of x. When nvar > 1, each of the
  % nvar equal-sized blacks of x are passed to func as a separate argument. 

  % Extract columns of the matrix into a cell array
  [n, ncols] = size (x);
  xcell = mat2cell (x, n, ncols / nvar * ones (1, nvar));

  % Evaluate column vectors as independent of arguments to bootfun
  retval = func (xcell{:});

end

%--------------------------------------------------------------------------

function [F, x] = empcdf (y, c)

  % Subfunction to calculate empirical cumulative distribution function
  %
  % Set c to:
  %  1 to have a complete distribution with F ranging from 0 to 1
  %  0 to avoid duplicate values in x
  %
  % Unlike ecdf, empcdf uses a denominator of N+1

  % Check input argument
  if ~isa(y,'numeric')
    error ('bootknife:empcdf: y must be numeric');
  end
  if all(size(y)>1)
    error ('bootknife:empcdf: y must be a vector');
  end
  if size(y,2)>1
    y = y.';
  end

  % Set defaults for optional input arguments
  if (nargin < 2)
    c = 1;
  end

  % Create empirical CDF
  ys = sort (y);
  N = sum (~ isnan (ys));
  [x, F] = unique (ys, 'last');
  F = F / (N + 1);

  % Apply option to complete the CDF
  if c > 0
    x = [x(1); x; x(end)];
    F = [0;F;1];
  end

  % Remove impossible values
  F(isnan (x)) = [];
  x(isnan (x)) = [];
  F(isinf (x)) = [];
  x(isinf (x)) = [];

end

%--------------------------------------------------------------------------

function X = kdeinv (P, Y, BW, CF)

  % Inverse of the cumulative density function (CDF) of a kernel density 
  % estimate (KDE)
  % 
  % The function returns X, the inverse CDF of the KDE of Y for the bandwidth
  % BW evaluated at the values in P. CF is a shrinkage factor for the variance
  % of the data in Y

  % Set defaults for optional input arguments
  if (nargin < 4)
    CF = 1;
  end

  % Create Normal CDF function
  pnorm = @(X, MU, SD) (0.5 * (1 + erf ((X - MU) / (SD * sqrt (2)))));

  % Calculate statistics of the data
  N = numel (Y);
  MU = mean (Y);
  SD = std (Y);

  % Apply shrinkage correction
  Y = ((Y - MU) * sqrt (CF)) + MU;

  % Set lower and upper bounds of X0 to range of Y
  X0 = [min(Y), max(Y)];

  % Perform root finding to get quantiles of the KDE at values of P
  findroot = @(X0, P) fzero (@(X) sum (pnorm (X - Y, 0, BW)) / N - P, X0);
  X = arrayfun (findroot, X0, P);

end

%--------------------------------------------------------------------------

function PX = ExpandProbs (P, DF)

  % Modify ALPHA to adjust tail probabilities assuming that the kurtosis of the
  % sampling distribution scales with degrees of freedom like the t-distribution.
  % This is related in concept to ExpandProbs in the resample R package:
  % https://www.rdocumentation.org/packages/resample/versions/0.6/topics/ExpandProbs

  % Get size of P
  sz = size (P);

  % Create required distribution functions
  stdnormcdf = @(X) 0.5 * (1 + erf (X / sqrt (2)));
  stdnorminv = @(P) sqrt (2) * erfinv (2 * P - 1);
  if exist('betaincinv','file')
    studinv = @(P, DF) sign (P - 0.5) * ...
                       sqrt ( DF ./ betaincinv (2 * min (P, 1 - P), DF / 2, 0.5) - DF);
  else
    % Earlier versions of matlab do not have betaincinv
    % Instead, use betainv from the Statistics and Machine Learning Toolbox
    try 
      studinv = @(P, DF) sign (P - 0.5) * ...
                         sqrt ( DF ./ betainv (2 * min (P, 1 - P), DF / 2, 0.5) - DF);
    catch
      % Use the Normal distribution (i.e. do not expand probabilities) if
      % either betaincinv or betainv are not available
      studinv = @(P, DF) stdnorminv (P);
      warning ('Could not create studinv function. Intervals will not be expanded.');
    end
  end
 
  % Calculate statistics of the data
  PX = stdnormcdf (arrayfun (studinv, P, repmat (DF, sz)));

end

%--------------------------------------------------------------------------

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! ## 95% BCa bootstrap confidence intervals for the mean
%! bootknife (data, 2000, @mean);

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! ## 95% calibrated percentile bootstrap confidence intervals for the mean
%! ## Calibration is of interval endpoints
%! bootknife (data, [2000, 200], @mean);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! ## 95% calibrated percentile bootstrap confidence intervals for the median
%! ## with smoothing. Calibration is of interval endpoints
%! bootknife (data, [2000, 200], @smoothmedian);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! ## 90% percentile bootstrap confidence intervals for the variance
%! bootknife (data, 2000, {@var, 1}, [0.05, 0.95]);

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! ## 90% BCa bootstrap confidence intervals for the variance
%! bootknife (data, 2000, {@var, 1}, 0.1);

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! ## 90% calibrated percentile bootstrap confidence intervals for the variance
%! ## Calibration of central coverage (2-sided)
%! bootknife (data, [2000, 200], {@var, 1}, 0.1);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! ## 90% calibrated percentile bootstrap confidence intervals for the variance
%! ## Calibration of central coverage (1-sided)
%! bootknife (data, [2000, 200], {@var, 1}, [0.05, 0.95]);
%!
%! ## Please be patient, the calculations will be completed soon...


%!demo
%!
%! ## Input bivariate dataset
%! x = [576 635 558 578 666 580 555 661 651 605 653 575 545 572 594]';
%! y = [3.39 3.3 2.81 3.03 3.44 3.07 3 3.43 3.36 3.13 3.12 2.74 2.76 2.88 2.96]'; 
%!
%! ## 95% BCa bootstrap confidence intervals for the correlation coefficient
%! bootknife ({x, y}, 2000, @corr);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Spatial Test Data from Table 14.1 of Efron and Tibshirani (1993)
%! ## An Introduction to the Bootstrap in Monographs on Statistics and Applied 
%! ## Probability 57 (Springer)
%!
%! ## AIM: to construct 90% nonparametric bootstrap confidence intervals for var(A,1)
%! ## var(A,1) = 171.5, and exact intervals based on Normal theory are [118.4, 305.2].
%! ## i.e. (numel(A)-1)*var(A,0) ./ chi2inv(1-[0.05;0.95],numel(A)-1)
%!
%! ## Calculations using the 'boot' and 'bootstrap' packages in R
%! ## 
%! ## library (boot)       # Functions from Davison and Hinkley (1997)
%! ## A <- c(48,36,20,29,42,42,20,42,22,41,45,14,6,0,33,28,34,4,32,24,47,41,24,26,30,41);
%! ## n <- length(A)
%! ## var.fun <- function (d, i) { 
%! ##               # Function to compute the population variance
%! ##               n <- length (d); 
%! ##               return (var (d[i]) * (n - 1) / n) };
%! ## boot.fun <- function (d, i) {
%! ##               # Compute the estimate
%! ##               t <- var.fun (d, i);
%! ##               # Compute sampling variance of the estimate using Tukey's jackknife
%! ##               n <- length (d);
%! ##               U <- empinf (data=d[i], statistic=var.fun, type="jack", stype="i");
%! ##               var.t <- sum (U^2 / (n * (n - 1)));
%! ##               return ( c(t, var.t) ) };
%! ## set.seed(1)
%! ## var.boot <- boot (data=A, statistic=boot.fun, R=20000, sim='balanced')
%! ## ci1 <- boot.ci (var.boot, conf=0.90, type="norm")
%! ## ci2 <- boot.ci (var.boot, conf=0.90, type="perc")
%! ## ci3 <- boot.ci (var.boot, conf=0.90, type="basic")
%! ## ci4 <- boot.ci (var.boot, conf=0.90, type="bca")
%! ## ci5 <- boot.ci (var.boot, conf=0.90, type="stud")
%! ##
%! ## library (bootstrap)  # Functions from Efron and Tibshirani (1993)
%! ## set.seed(1); ci4a <- bcanon (A, 20000, var.fun, alpha=c(0.05,0.95))
%! ## set.seed(1); ci5a <- boott (A, var.fun, nboott=20000, nbootsd=500, perc=c(.05,.95))
%! ##
%! ## Summary of confidence intervals from 'boot' and 'bootstrap' packages in R
%! ##
%! ## method                      |   0.05 |   0.95 | length | shape |  
%! ## ----------------------------|--------|--------|--------|-------|
%! ## ci1  - normal               |  109.4 |  246.8 |  137.4 |  1.21 |
%! ## ci2  - percentile           |   97.8 |  235.6 |  137.8 |  0.87 |
%! ## ci3  - basic                |  107.4 |  245.3 |  137.9 |  1.15 |
%! ## ci4  - BCa                  |  116.9 |  264.0 |  147.1 |  1.69 |
%! ## ci4a - BCa                  |  116.2 |  264.0 |  147.8 |  1.67 |
%! ## ci5  - bootstrap-t          |  111.8 |  291.2 |  179.4 |  2.01 |
%! ## ci5a - bootstrap-t          |  112.7 |  292.6 |  179.9 |  2.06 |
%! ## ----------------------------|--------|--------|--------|-------|
%! ## parametric - exact          |  118.4 |  305.2 |  186.8 |  2.52 |
%! ##
%! ## Summary of bias statistics from 'boot' package in R
%! ##
%! ## method                      | original |    bias | bias-corrected |
%! ## ----------------------------|----------|---------|----------------|
%! ## single bootstrap            |   171.53 |   -6.62 |         178.16 |
%! ## ----------------------------|----------|---------|----------------|
%! ## parametric - exact          |   171.53 |   -6.86 |         178.40 |
%! 
%! ## Calculations using the 'statistics-bootstrap' package for Octave/Matlab
%! ##
%! ## A = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%! ##      0 33 28 34 4 32 24 47 41 24 26 30 41]';
%! ## boot (1,1,false,1); ci2 = bootknife (A, 20000, {@var,1}, [0.05,0.95]);
%! ## boot (1,1,false,1); ci4 = bootknife (A, 20000, {@var,1}, 0.1);
%! ## boot (1,1,false,1); ci6a = bootknife (A, [20000,500], {@var,1}, 0.1);
%! ## boot (1,1,false,1); ci6b = bootknife (A, [20000,500], {@var,1}, [0.05,0.95]);
%! ##
%! ## Summary of confidence intervals from 'statistics-bootstrap' package for Octave/Matlab
%! ##
%! ## method                      |   0.05 |   0.95 | length | shape |
%! ## ----------------------------|--------|--------|--------|-------|
%! ## ci2  - percentile           |   96.3 |  237.0 |  140.7 |  0.87 |
%! ## ci4  - BCa                  |  115.3 |  263.0 |  147.7 |  1.63 |
%! ## ci6a - calibrated (2-sided) |   81.8 |  254.4 |  172.6 |  0.92 |
%! ## ci6b - calibrated (1-sided) |  114.3 |  295.5 |  181.2 |  2.17 |
%! ## ----------------------------|--------|--------|--------|-------|
%! ## parametric - exact          |  118.4 |  305.2 |  186.8 |  2.52 |
%! ##
%! ## Simulation results for constructing 90% confidence intervals for the
%! ## variance of a population N(0,1) from 1000 random samples of size 26
%! ## (analagous to the situation above). Simulation performed using the
%! ## bootsim script with nboot of 2000 (for single bootstrap) or [2000,200]
%! ## (for double bootstrap).
%! ##
%! ## method               | coverage |  lower |  upper | length | shape |
%! ## ---------------------|----------|--------|--------|--------|-------|
%! ## percentile           |    81.8% |   1.3% |  16.9% |   0.80 |  0.92 |
%! ## BCa                  |    87.3% |   4.2% |   8.5% |   0.86 |  1.85 |
%! ## calibrated (2-sided) |    90.5% |   0.5% |   9.0% |   1.06 |  1.06 |
%! ## calibrated (1-sided) |    90.7% |   5.1% |   4.2% |   1.13 |  2.73 |
%! ## ---------------------|----------|--------|--------|--------|-------|
%! ## parametric - exact   |    90.8% |   3.7% |   5.5% |   0.99 |  2.52 |
%!
%! ## Summary of bias statistics from 'boot' package in R
%! ##
%! ## method             | original |    bias | bias-corrected |
%! ## -------------------|----------|---------|----------------|
%! ## single bootstrap   |   171.53 |   -6.70 |         178.24 |
%! ## double bootstrap   |   171.53 |   -6.85 |         178.38 |
%! ## -------------------|----------|---------|----------------|
%! ## parametric - exact |   171.53 |   -6.86 |         178.40 |
%!
%! ## The equivalent methods for constructing bootstrap intervals in the 'boot'
%! ## and 'bootstrap' packages (in R) and the statistics-bootstrap package (in
%! ## Octave/Matlab) produce intervals with very similar end points, length and
%! ## shape. However, all intervals calculated using the 'statistics-bootstrap'
%! ## package are slightly longer than the intervals calculated in R because
%! ## the 'statistics-bootstrap' package uses bootknife resampling. The scale of
%! ## the sampling distribution for small samples is approximated better by
%! ## bootknife (rather than bootstrap) resampling. 

%!test
%! ## Spatial test data from Table 14.1 of Efron and Tibshirani (1993)
%! ## An Introduction to the Bootstrap in Monographs on Statistics and Applied 
%! ## Probability 57 (Springer)
%! A = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!      0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! ## Nonparametric 90% percentile confidence intervals
%! ## Table 14.2 percentile intervals are 100.8 - 233.9
%! boot (1, 1, false, 1); # Set random seed
%! stats = bootknife(A,2000,{@var,1},[0.05 0.95]);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (stats.original, 171.534023668639, 1e-09);
%!   assert (stats.bias, -7.323387573964482, 1e-09);
%!   assert (stats.std_error, 43.30079972388541, 1e-09);
%!   assert (stats.CI_lower, 95.24158837039771, 1e-09);
%!   assert (stats.CI_upper, 237.7156378257705, 1e-09);
%! end
%!
%! ## Nonparametric 90% BCa confidence intervals
%! ## Table 14.2 BCa intervals are 115.8 - 259.6
%! boot (1, 1, false, 1); # Set random seed
%! stats = bootknife(A,2000,{@var,1},0.1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (stats.original, 171.534023668639, 1e-09);
%!   assert (stats.bias, -7.323387573964482, 1e-09);
%!   assert (stats.std_error, 43.30079972388541, 1e-09);
%!   assert (stats.CI_lower, 113.2388308884533, 1e-09);
%!   assert (stats.CI_upper, 264.9901439787903, 1e-09);
%! end
%!
%! ## Nonparametric 90% calibrated percentile confidence intervals (2-sided)
%! boot (1, 1, false, 1); # Set random seed
%! stats = bootknife(A,[2000,200],{@var,1},0.1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (stats.original, 171.534023668639, 1e-09);
%!   assert (stats.bias, -8.088193809171344, 1e-09);
%!   assert (stats.std_error, 46.53418481731099, 1e-09);
%!   assert (stats.CI_lower, 78.7080103220452, 1e-09);
%!   assert (stats.CI_upper, 260.8289943527518, 1e-09);
%! end
%!
%! ## Nonparametric 90% calibrated percentile confidence intervals (1-sided)
%! boot (1, 1, false, 1); # Set random seed
%! stats = bootknife(A,[2000,200],{@var,1},[0.05,0.95]);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (stats.original, 171.534023668639, 1e-09);
%!   assert (stats.bias, -8.088193809171344, 1e-09);
%!   assert (stats.std_error, 46.53418481731099, 1e-09);
%!   assert (stats.CI_lower, 110.7774121402538, 1e-09);
%!   assert (stats.CI_upper, 303.7824606292835, 1e-09);
%! end
%!
%! ## Exact intervals based on normal theory are 118.4 - 305.2 (Table 14.2)
%! ## Note that all of the bootknife intervals are slightly wider than the
%! ## non-parametric intervals in Table 14.2 because the bootknife (rather than
%! ## standard bootstrap) resampling used here reduces small sample bias

%!test
%! ## Law school data from Table 3.1 of Efron and Tibshirani (1993)
%! ## An Introduction to the Bootstrap in Monographs on Statistics and Applied 
%! ## Probability 57 (Springer)
%! LSAT = [576 635 558 578 666 580 555 661 651 605 653 575 545 572 594]';
%! GPA = [3.39 3.3 2.81 3.03 3.44 3.07 3 3.43 3.36 3.13 3.12 2.74 2.76 2.88 2.96]';
%!
%! ## Nonparametric 90% percentile confidence intervals
%! ## Percentile intervals on page 266 are 0.524 - 0.928
%! boot (1, 1, false, 1); # Set random seed
%! stats = bootknife({LSAT,GPA},2000,@cor,[0.05,0.95]);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (stats.original, 0.7763744912894069, 1e-09);
%!   assert (stats.bias, -0.008259337758776852, 1e-09);
%!   assert (stats.std_error, 0.1420949476115542, 1e-09);
%!   assert (stats.CI_lower, 0.5056363801008388, 1e-09);
%!   assert (stats.CI_upper, 0.9586254199016858, 1e-09);
%! end
%!
%! ## Nonparametric 90% BCa confidence intervals
%! ## BCa intervals on page 266 are 0.410 - 0.923
%! boot (1, 1, false, 1); # Set random seed
%! stats = bootknife({LSAT,GPA},2000,@cor,0.1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (stats.original, 0.7763744912894069, 1e-09);
%!   assert (stats.bias, -0.008259337758776852, 1e-09);
%!   assert (stats.std_error, 0.1420949476115542, 1e-09);
%!   assert (stats.CI_lower, 0.4119228032301614, 1e-09);
%!   assert (stats.CI_upper, 0.9300646701004258, 1e-09);
%! end
%!
%! ## Nonparametric 90% calibrated percentile confidence intervals (1-sided)
%! boot (1, 1, false, 1); # Set random seed
%! stats = bootknife({LSAT,GPA},[2000,500],@cor,[0.05,0.95]);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (stats.original, 0.7763744912894069, 1e-09);
%!   assert (stats.bias, -0.009420108365347568, 1e-09);
%!   assert (stats.std_error, 0.1438249935781226, 1e-09);
%!   assert (stats.CI_lower, 0.254445053627807, 1e-09);
%!   assert (stats.CI_upper, 0.9469388762553252, 1e-09);
%! end
%! ## Exact intervals based on normal theory are 0.51 - 0.91