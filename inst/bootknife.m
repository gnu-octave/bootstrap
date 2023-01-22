%  Function File: bootknife
%
%  Bootknife (bootstrap) resampling and statistics
%
%  This function takes a DATA sample (containing n rows) and uses bootstrap
%  resampling methods to calculate a bias of the parameter estimate, a standard
%  error, and confidence intervals. Specifically, the method uses bootknife
%  resampling [1], which involves creating leave-one-out jackknife samples
%  of size n - 1 and then drawing samples of size n with replacement from the
%  jackknife samples. The resampling of DATA rows is balanced in order to
%  reduce Monte Carlo error [2,3].
%
%  If single bootstrap is requested, the confidence intervals are obtained from 
%  the quantiles of a kernel density estimate of the bootstrap statistics (with
%  shrinkage corrrection). By default, the confidence intervals are bias-
%  corrected and accelerated (BCa) [4-5]. BCa intervals are fast to compute and
%  have good coverage and correctness when combined with bootknife resampling
%  as it is here [1].
%
%  If double bootstrap is requested, the algorithm uses calibration to improve
%  the accuracy (of the bias and standard error) and coverage of the confidence
%  intervals [6-12], which are obtained from the empirical distribution of the
%  bootstrap statistics by linear interpolation [9]. 
%
%  STATS = bootknife (DATA)
%  STATS = bootknife ({DATA})
%  STATS = bootknife (DATA, NBOOT)
%  STATS = bootknife (DATA, NBOOT, BOOTFUN)
%  STATS = bootknife (DATA, NBOOT, {BOOTFUN, ...})
%  STATS = bootknife (DATA, NBOOT, ..., ALPHA)
%  STATS = bootknife (DATA, NBOOT, ..., ALPHA, STRATA)
%  STATS = bootknife (DATA, NBOOT, ..., ALPHA, STRATA, NPROC)
%  STATS = bootknife (DATA, NBOOT, ..., ALPHA, STRATA, NPROC, BOOTSAM)
%  [STATS, BOOTSTAT] = bootknife (...)
%  [STATS, BOOTSTAT] = bootknife (...)
%  [STATS, BOOTSTAT, BOOTSAM] = bootknife (...)
%  bootknife (DATA,...);
%  bootknife (DATA, [2000, 0], @mean, 0.05, [], 0)            % Defaults (single)
%  bootknife (DATA, [2000, 200], @mean, [0.025,0.975], [], 0) % Defaults (double)
%
%  STATS = bootknife (DATA) resamples from the rows of a DATA sample (column 
%  vector or matrix) and returns a structure (or structure array) with the
%  following fields:
%    original: contains the result of applying BOOTFUN to the DATA 
%    bias: contains the bootstrap estimate of bias [11-12]
%    std_error: contains the bootstrap estimate of the standard error of bootfun
%    CI_lower: contains the lower bound of the bootstrap confidence interval
%    CI_upper: contains the upper bound of the bootstrap confidence interval
%  By default, the statistics relate to BOOTFUN being @mean and the confidence
%  intervals are 95% bias-corrected and accelerated (BCa) intervals [1,4-5]. 
%
%  STATS = bootknife ({DATA}) resamples from the rows of DATA as above. If
%  DATA (in the cell array) is multiple column vectors and/or matrices, then
%  bootknife resamples the rows across all the column vectors or matrices, and
%  the resamples are passed to BOOTFUN as multiple data input arguments. All
%  DATA column vectors and matrices must have the same number of rows.
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
%  have nominal probabilities of 2.5% and 97.5% [5]. Note that one can get
%  away with a lower number of resamples in the second bootstrap to reduce
%  the computational expense of the double bootstrap (e.g. [2000, 200]),
%  since the algorithm uses linear interpolation to achieve near-asymptotic
%  calibration of confidence intervals [3]. The confidence intervals calculated
%  (with either single or double bootstrap) are transformation respecting and
%  can have better coverage and/or accuracy compared to intervals derived from
%  normal theory or to simple percentile bootstrap confidence intervals [5-10].
%
%  STATS = bootknife (DATA, NBOOT, BOOTFUN) also specifies BOOTFUN, a function 
%  handle, a string indicating the name of the function to apply to the DATA
%  (and each bootstrap resample), or a cell array where the first cell is the 
%  function handle or string, and other cells being additional input arguments 
%  for BOOTFUN, where BOOTFUN must take DATA for the first input arguments.
%  BOOTFUN can return a scalar or any multidimensional numeric variable, but
%  the output will be reshaped as a column vector. BOOTFUN must calculate a
%  statistic representative of the finite DATA sample, it should not be an
%  unbiased estimate of a population parameter. For example, for the variance,
%  set BOOTFUN to {@var,1} (not @var or {@var,0}), and for the correlation
%  coefficient, use the provided @cor function. Smooth functions of the DATA
%  are preferable, (e.g. use smoothmedian function instead of the ordinary
%  median). The default value(s) of BOOTFUN is/are the (column) mean(s).
%    When BOOTFUN is @mean or 'mean', residual narrowness bias of central
%  coverage is almost eliminated by using Student's t-distribution to expand
%  the probabilities of percentiles [13].
%    Note that if single bootstrap is requested and BOOTFUN cannot be executed
%  during leave-one-out jackknife, the acceleration constant will be set to 0
%  and intervals will be bias-corrected only.
%    Recommended BOOFUN for some commonly used statistics:
%      - Mean: @mean
%      - Standard deviation: {@std,1}
%      - Variance: {@var,1}
%      - Correlation coefficient: @cor
%      - Linear regression: @(y,X) pinv(X)*y or @regress
%      - Median: @median or @smoothmedian
%      - Median absolute deviation: @mad or @smoothmad
%    See code comments or demos for examples of usage.
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
%  - PERCENTILE (equal-tailed): ALPHA must be a pair of probabilities and NBOOT
%    must be a scalar value (or the second element in NBOOT must be zero).
%
%  - BIAS-CORRECTED AND ACCELERATED (BCa): ALPHA must be a scalar value and
%    NBOOT must be a scalar value (or the second element in NBOOT must be zero).
%
%  - CALIBRATED PERCENTILE (equal-tailed): ALPHA must be a scalar value and
%    NBOOT must be a vector of two positive, non-zero integers (for double
%    bootstrap). The interval construction corresponds to the 2-sided intervals
%    intervals in [7,9].
%
%  - CALIBRATED PERCENTILE: ALPHA must be must be a pair of probabilities and
%    NBOOT must be a vector of two positive, non-zero integers (for double
%    bootstrap). The interval construction corresponds to the lower and upper
%    1-sided intervals in [7-10]. The method used is equivalent to the
%    confidence point calibration algorithm 18.1 in [5].
%
%  Confidence intervals are not calculated when the value(s) of ALPHA is/are
%  NaN. If empty (or not specified), the default value of ALPHA is 0.05 for
%  single bootstrap (i.e. 95% BCa intervals), or [0.025, 0.975] for double 
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
%  STATS = bootknife (DATA, NBOOT, ..., ALPHA, STRATA, NPROC, BOOTSAM) uses
%  bootstrap resampling indices provided in BOOTSAM. The BOOTSAM should be a
%  matrix with the same number of rows as the data. When BOOTSAM is provided,
%  the first element of NBOOT is ignored.
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
%  [7] Lee and Young (1999) The effect of Monte Carlo approximation on coverage
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
%  bootknife (version 2023.01.19)
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


function [stats, bootstat, bootsam] = bootknife (x, nboot, bootfun, alpha, ...
                              strata, ncpus, bootsam, REF, ISOCTAVE, ERRCHK)

  % Input argument names in all-caps are for internal use only
  % REF, ISOCTAVE and ERRCHK are undocumented input arguments required
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
    if ((nargin < 2) || isempty (nboot))
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
      bootfun_str = 'mean';
    else
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
        bootfun = str2func (bootfun);
      elseif (isa (bootfun, 'function_handle'))
        bootfun_str = func2str (bootfun);
      else
        error ('bootknife: BOOTFUN must be a function name or function handle')
      end
    end
    if (iscell (x))
      % If DATA is a cell array of equal size colunmn vectors, convert the cell
      % array to a matrix and redefine bootfun to parse multiple input arguments
      szx = cellfun (@(x) size (x,2), x);
      x = [x{:}];
      bootfun = @(x) localfunc.col2args (bootfun, x, szx);
    else
      szx = size (x, 2);
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
      if (~ isa (alpha, 'numeric') || (nalpha > 2))
        error ('bootknife: ALPHA must be a scalar (two-tailed probability) or a vector (pair of probabilities)');
      end
      if (size (alpha, 1) > 1)
        alpha = alpha.';
      end
      if (any ((alpha < 0) | (alpha > 1)))
        error ('bootknife: Value(s) in ALPHA must be between 0 and 1');
      end
      if (nalpha > 1)
        % alpha is a pair of probabilities
        % Make sure probabilities are in the correct order
        if (alpha(1) > alpha(2) )
          error ('bootknife: The pair of probabilities must be in ascending numeric order');
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
      if (~ isa (ncpus, 'numeric'))
        error ('bootknife: NPROC must be numeric');
      end
      if (any (ncpus ~= abs (fix (ncpus))))
        error ('bootknife: NPROC must be a positive integer');
      end    
      if (numel (ncpus) > 1)
        error ('bootknife: NPROC must be a scalar value');
      end
    end
    if ((nargin < 9) || isempty (ISOCTAVE))
      % Check if running in Octave (else assume Matlab)
      info = ver; 
      ISOCTAVE = any (ismember ({info.Name}, 'Octave'));
    end
    if (ISOCTAVE)
      ncpus = min (ncpus, nproc);
    else
      ncpus = min (ncpus, feature ('numcores'));
    end
  else
    szx = 1;
  end

  % Determine properties of the DATA (x)
  [n, nvar] = size (x);
  if (n < 2)
    error ('bootknife: DATA must be numeric and contain > 1 row')
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
  if (any (isnan (T0)))
    error ('bootknife: BOOTFUN returned NaN with the DATA provided')
  end

  % Check whether bootfun is vectorized
  if (nvar > 1)
    M = cell2mat (cellfun (@(i) repmat (x(:, i), 1, 2), ...
                  num2cell (1:nvar), 'UniformOutput', false));
  else
    M = repmat (x, 1, 2);
  end
  if (any (szx > 1))
    vectorized = false;
  else
    try
      chk = bootfun (M);
      if (all (size (chk) == [size(T0, 1), 2]) && all (chk == bootfun (x)))
        vectorized = true;
      else
        vectorized = false;
      end
    catch
      vectorized = false;
    end
  end

  % Initialize probabilities
  l = [];

  % If applicable, check we have parallel computing capabilities
  if (ncpus > 1)
    if (ISOCTAVE)
      pat = '^parallel';
      software = pkg ('list');
      names = cellfun (@(S) S.name, software, 'UniformOutput', false);
      status = cellfun (@(S) S.loaded, software, 'UniformOutput', false);
      index = find (~ cellfun (@isempty, regexpi (names,pat)));
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
      info = ver; 
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
        warning ('bootknife:parallel', ...
           'Parallel Computing Toolbox not installed or operational. Falling back to serial processing.');
        ncpus = 1;
      end
    end
  else
    if ((ncpus > 1) && ~ PARALLEL)
      if (ISOCTAVE)
        % OCTAVE Parallel Computing Package is not installed or loaded
        warning ('bootknife:parallel', ...
          'Parallel Computing Package not installed and/or loaded. Falling back to serial processing.');
      else
        % MATLAB Parallel Computing Toolbox is not installed or loaded
        warning ('bootknife:parallel', ...
          'Parallel Computing Toolbox not installed and/or loaded. Falling back to serial processing.');
      end
      ncpus = 0;
    end
  end

  % Calculate the number of elements in the return value of bootfun 
  m = numel (T0);
  if (m > 1)
    % Vectorized along the dimension of the return values of bootfun so
    % reshape the output to be a column vector before proceeding with bootstrap
    if (size (T0, 2) > 1)
      bootfun = @(x) reshape (bootfun (x), [], 1);
      T0 = reshape (T0, [], 1);
      vectorized = false;
    end
  end

  % Evaluate strata input argument
  if (~ isempty (strata))
    if (~ isnumeric (strata))
      % Convert strata to numeric ID
      [jnk1, jnk2, strata] = unique (strata);
      clear jnk1 jnk2;
    end
    % Get strata IDs
    gid = unique (strata);  % strata ID
    K = numel (gid);        % number of strata
    % Create strata matrix
    g = false (n, K);
    for k = 1:K
      g(:, k) = (strata == gid(k));
    end
    nk = sum (g);          % strata sample sizes
  else 
    g = ones (n, 1);
    K = 1;
  end

  % Perform balanced bootknife resampling
  unbiased = true;  % Set to true for bootknife resampling
  if ((nargin < 7) || isempty (bootsam))
    if (~ isempty (strata))
      if (nvar > 1) || (nargout > 2)
        % We can save some memory by making bootsam an int32 datatype
        bootsam = zeros (n, B, 'int32');
        for k = 1:K
          if ((sum (g(:, k))) > 1)
            bootsam(g(:, k), :) = boot (find (g(:, k)), B, unbiased);
          else
            bootsam(g(:, k), :) = find (g(:, k)) * ones (1, B);
          end
        end
      else
        % For more efficiency, if we don't need bootsam, we can directly resample values of x
        bootsam = [];
        X = zeros (n, B);
        for k = 1:K
          if ((sum (g(:, k))) > 1)
            X(g(:, k), :) = boot (x(g(:, k), :), B, unbiased);
          else
            X(g(:, k), :) = x(g(:, k), :) * ones (1, B);
          end
        end
      end
    else
      if (nvar > 1) || (nargout > 2)
        % We can save some memory by making bootsam an int32 datatype
        bootsam = zeros (n, B, 'int32');
        bootsam(:, :) = boot (n, B, unbiased);
      else
        % For more efficiency, if we don't need bootsam, we can directly resample values of x
        bootsam = [];
        X = boot (x, B, unbiased);
      end
    end
  else
    if (size (bootsam, 1) ~= n)
      error ('bootknife: BOOTSAM must have the same number of rows as X')
    end
    nboot(1) = size (bootsam, 2);
    B = nboot(1);
  end

  % Evaluate bootfun each bootstrap resample
  if (isempty (bootsam))
    if (vectorized)
      % Vectorized evaluation of bootfun on the DATA resamples
      bootstat = bootfun (X);
    else
      if (ncpus > 1)
        % Evaluate bootfun on each bootstrap resample in PARALLEL
        if (ISOCTAVE)
          % OCTAVE
          bootstat = parcellfun (ncpus, bootfun, num2cell (X, 1), 'UniformOutput', false);
        else
          % MATLAB
          bootstat = cell (1, B);
          parfor b = 1:B; bootstat{b} = bootfun (X(:, b)); end
        end
      else
        bootstat = cellfun (bootfun, num2cell (X, 1), 'UniformOutput', false);
      end
    end
  else
    if (vectorized)
      % DATA resampling (using bootsam) and vectorized evaluation of bootfun on 
      % the DATA resamples 
      if (nvar > 1)
        % Multivariate
        % Perform DATA sampling
        X = cell2mat (cellfun (@(i) reshape (x(bootsam, i), n, B), ...
                      num2cell (1:nvar, 1), 'UniformOutput', false));
      else
        % Univariate
        % Perform DATA sampling
        X = x(bootsam);
      end
      % Function evaluation on bootknife samples
      bootstat = bootfun (X);
    else 
      cellfunc = @(bootsam) bootfun (x(bootsam, :));
      if (ncpus > 1)
        % Evaluate bootfun on each bootstrap resample in PARALLEL
        if (ISOCTAVE)
          % OCTAVE
          bootstat = parcellfun (ncpus, cellfunc, num2cell (bootsam, 1), 'UniformOutput', false);
        else
          % MATLAB
          bootstat = cell (1, B);
          parfor b = 1:B; bootstat{b} = cellfunc (bootsam(:, b)); end
        end
      else
        % Evaluate bootfun on each bootstrap resample in SERIAL
        bootstat = cellfun (cellfunc, num2cell (bootsam, 1), 'UniformOutput', false);
      end
    end
  end
  if (iscell (bootstat))
    bootstat = cell2mat (bootstat);
  end
  
  % Remove bootstrap statistics that contain NaN, along with their associated 
  % DATA resamples in X or bootsam
  ridx = any (isnan (bootstat), 1);
  bootstat_all = bootstat;
  bootstat(:, ridx) = [];
  if (isempty (bootsam))
    X(:, ridx) = [];
  else
    bootsam(:, ridx) = [];
  end
  if (isempty (bootstat))
    error ('bootknife: BOOTFUN returned NaN for every bootstrap resamples')
  end
  B = B - sum (ridx);

  % Calculate the bootstrap bias, standard error and confidence intervals 
  if (C > 0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%% DOUBLE BOOTSTRAP %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (ncpus > 1)
      % PARALLEL execution of inner layer resampling for double (i.e. iterated) bootstrap
      if (ISOCTAVE)
        % OCTAVE
        % Set unique random seed for each parallel thread
        pararrayfun (ncpus, @boot, 1, 1, false, 1:ncpus);
        if (vectorized && isempty (bootsam))
          cellfunc = @(x) bootknife (x, C, bootfun, NaN, strata, 0, [], T0, ISOCTAVE, false);
          bootout = parcellfun (ncpus, cellfunc, num2cell (X, 1), 'UniformOutput', false);
        else
          cellfunc = @(bootsam) bootknife (x(bootsam, :), C, bootfun, NaN, strata, 0, [], T0, ISOCTAVE, false);
          bootout = parcellfun (ncpus, cellfunc, num2cell (bootsam, 1), 'UniformOutput', false);
        end
      else
        % MATLAB
        % Set unique random seed for each parallel thread
        parfor i = 1:ncpus; boot (1, 1, false, i); end
        % Perform inner layer of resampling
        % Preallocate structure array
        bootout = cell (1, B);
        if (vectorized && isempty (bootsam))
          cellfunc = @(x) bootknife (x, C, bootfun, NaN, strata, 0, [], T0, ISOCTAVE, false);
          parfor b = 1:B; bootout{b} = cellfunc (X(:, b)); end
        else
          cellfunc = @(bootsam) bootknife (x(bootsam, :), C, bootfun, NaN, strata, 0, [], T0, ISOCTAVE, false);
          parfor b = 1:B; bootout{b} = cellfunc (bootsam(:, b)); end
        end
      end
    else
      % SERIAL execution of inner layer resampling for double bootstrap
      if (vectorized && isempty (bootsam))
        cellfunc = @(x) bootknife (x, C, bootfun, NaN, strata, 0, [], T0, ISOCTAVE, false);
        bootout = cellfun (cellfunc, num2cell (X, 1), 'UniformOutput', false);
      else
        cellfunc = @(bootsam) bootknife (x(bootsam, :), C, bootfun, NaN, strata, 0, [], T0, ISOCTAVE, false);
        bootout = cellfun (cellfunc, num2cell (bootsam, 1), 'UniformOutput', false);
      end
    end
    % Double bootstrap bias estimation
    mu = cell2mat (cellfun (@(S) S.bias, bootout, 'UniformOutput', false)) + ...
         cell2mat (cellfun (@(S) S.original, bootout, 'UniformOutput', false));
    b = mean (bootstat, 2) - T0;
    c = mean (mu, 2) - 2 * mean (bootstat, 2) + T0;
    bias = b - c;
    % Double bootstrap multiplicative correction of the standard error
    V = cell2mat (cellfun (@(S) S.std_error.^2, bootout, 'UniformOutput', false));
    se = sqrt (var (bootstat, 0, 2).^2 ./ mean (V, 2));
    % Double bootstrap confidence intervals
    if (~ isnan (alpha))
      U = cell2mat (cellfun (@(S) S.Pr, bootout, 'UniformOutput', false));
      l = zeros (m, 2);
      ci = zeros (m, 2);
      for j = 1:m
        % Calibrate interval coverage
        switch nalpha
          case 1
            % alpha is a two-tailed probability (scalar)
            % Calibrate central coverage and construct equal-tailed intervals (2-sided)
            [cdf, v] = localfunc.empcdf (abs (2 * U(j, :) - 1));
            vk = interp1 (cdf, v, 1 - alpha, 'linear');
            l(j, :) = arrayfun (@(sign) 0.5 * (1 + sign * vk), [-1, 1]);
          case 2
            % alpha is a pair of probabilities (vector)
            % Calibrate coverage but construct endpoints separately (1-sided)
            % This is equivalent to algorithm 18.1 in Efron, and Tibshirani (1993)
            [cdf, u] = localfunc.empcdf (U(j, :));
            l(j, :) = arrayfun (@(p) interp1 (cdf, u, p, 'linear'), alpha);
        end
        % Linear interpolation
        [cdf, t1] = localfunc.empcdf (bootstat(j, :));
        ci(j, :) = arrayfun (@(p) interp1 (cdf, t1, p, 'linear'), l(j, :));
      end
    else
      ci = nan (m, 2);
    end
  else
    %%%%%%%%%%%%%%%%%%%%%%%%%%% SINGLE BOOTSTRAP %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bootstrap bias estimation
    bias = mean (bootstat, 2) - T0;
    % Bootstrap standard error
    se = std (bootstat, 0, 2);  % Unbiased estimate since we used bootknife resampling
    if (~ isnan (alpha))
      % If bootfun is the arithmetic meam, perform expansion based on t-distribution
      if (strcmpi (bootfun_str, 'mean'))
        expan_alpha = (3 - nalpha) * localfunc.ExpandProbs (alpha / (3 - nalpha), n - K);
      else
        expan_alpha = alpha;
      end
      state = warning;
      if (ISOCTAVE)
        warning ('on', 'quiet');
      else
        warning ('off', 'all');
      end
      % Create distribution functions
      stdnormcdf = @(x) 0.5 * (1 + erf (x / sqrt (2)));
      stdnorminv = @(p) sqrt (2) * erfinv (2 * p-1);
      switch nalpha
        case 1
          % Attempt to form bias-corrected and accelerated (BCa) bootstrap confidence intervals. 
          % Use the Jackknife to calculate the acceleration constant (a)
          try
            jackfun = @(i) bootfun (x(1:n ~= i, :));
            if (ncpus > 1)  
              % PARALLEL evaluation of bootfun on each jackknife resample 
              if (ISOCTAVE)
                % OCTAVE
                T = cell2mat (pararrayfun (ncpus, jackfun, 1:n, 'UniformOutput', false));
              else
                % MATLAB
                T = zeros (m, n);
                parfor i = 1:n; T(:, i) = feval (jackfun, i); end
              end
            else
              % SERIAL evaluation of bootfun on each jackknife resample
              T = cell2mat (arrayfun (jackfun, 1:n, 'UniformOutput', false));
            end
            % Calculate empirical influence function
            if (~ isempty (strata))
              gk = sum (g .* repmat (nk, n, 1), 2).';
              U = bsxfun (@times, gk - 1, bsxfun (@minus, mean (T, 2), T));  
            else
              U = (n - 1) * bsxfun (@minus, mean (T, 2), T);
            end
            a = sum (U.^3, 2) ./ (6 * sum (U.^2, 2) .^ 1.5);
          catch
            % Revert to bias-corrected (BC) bootstrap confidence intervals
            warning ('bootknife:jackfail', ...
              'BOOTFUN failed during jackknife calculations; acceleration constant set to 0\n');
            a = zeros (m, 1);
          end
          % Calculate the bias correction constant (z0)
          % Calculate the median bias correction z0
          z0 = stdnorminv (sum (bsxfun (@lt, bootstat, T0), 2) / B);
          if (~ all (isfinite (z0)))
            % Revert to percentile bootstrap confidence intervals
            warning ('bootknife:biasfail', ...
              'Unable to calculate the bias correction constant; reverting to percentile intervals\n');
            z0 = zeros (m, 1);
            a = zeros (m, 1); 
            l = repmat ([expan_alpha / 2, 1 - expan_alpha / 2], m, 1);
          end
          if (isempty (l))
            % Calculate BCa or BC percentiles
            z1 = stdnorminv (expan_alpha / 2);
            z2 = stdnorminv (1 - expan_alpha / 2);
            l = cat (2, stdnormcdf (z0 + ((z0 + z1) ./ (1 - bsxfun (@times, a , z0 + z1)))),... 
                        stdnormcdf (z0 + ((z0 + z2) ./ (1 - bsxfun (@times, a , z0 + z2)))));
          end
        case 2
          % alpha is a vector of probabilities
          l = repmat (expan_alpha, m, 1);
      end
      % Intervals constructed from kernel density estimate of the bootstrap
      % (with shrinkage correction)
      ci = zeros (m, 2);
      for j = 1:m
        try
          ci(j, :) = localfunc.kdeinv (l(j, :), bootstat(j, :), se(j) * sqrt (1 / (n - K)), 1 - 1 / (n - K));
        catch
          % Linear interpolation (legacy)
          fprintf ('Note: Falling back to linear interpolation to calculate percentiles for interval pair %u\n', j);
          [cdf, t1] = localfunc.empcdf (bootstat(j, :));
          ci(j, :) = arrayfun (@(p) interp1 (cdf, t1, p, 'linear'), l(j, :));
        end
      end
      warning (state);
      if (ISOCTAVE)
        warning ('off', 'quiet');
      end
    else
      ci = nan (m, 2);
    end
  end

  % Prepare output arguments
  stats = struct;
  stats.original = T0;
  stats.bias = bias;
  stats.std_error = se;
  stats.CI_lower = ci(:, 1);
  stats.CI_upper = ci(:, 2);
  % Use quick interpolation to find the proportion (Pr) of bootstat <= REF
  if ((nargin > 7) && ~ isempty (REF))
    I = bsxfun (@le, bootstat, REF);
    pr = sum (I, 2);
    t = cell2mat (arrayfun (@(j) ...
         [max([min(bootstat(j, :)), max(bootstat(j, I(j, :)))]),...
          min([max(bootstat(j, :)), min(bootstat(j, ~ I(j, :)))])], (1:m).', ...
          'UniformOutput', false));
    dt = t(:, 2) - t(:, 1);
    chk = and (pr < B, dt > 0);
    Pr = zeros (m, 1);
    Pr(chk, 1) = pr(chk, 1) + ((REF(chk, 1) - t(chk, 1) ).* ...
                        (min (pr(chk, 1) + 1, B) - pr(chk, 1)) ./ dt(chk, 1));
    Pr(~ chk, 1) = pr(~ chk, 1);
    stats.Pr = Pr / B;
  end
  bootstat = bootstat_all;

  % Print output if no output arguments are requested
  if (nargout == 0) 
    print_output (stats, nboot, alpha, l, m, bootfun_str, strata);
  else
    if (isempty (bootsam))
      [warnmsg, warnID] = lastwarn;
      if (ismember (warnID, {'bootknife:biasfail','bootknife:jackfail'}))
        warning ('bootknife:lastwarn', warnmsg);
      end
      lastwarn ('', '');
    end
  end

end

%--------------------------------------------------------------------------

function print_output (stats, nboot, alpha, l, m, bootfun_str, strata)

    fprintf (['\nSummary of non-parametric bootstrap estimates of bias and precision\n',...
              '******************************************************************************\n\n']);
    fprintf ('Bootstrap settings: \n');
    fprintf (' Function: %s\n', bootfun_str);
    if (nboot(2) > 0)
      if (isempty (strata))
        fprintf (' Resampling method: Iterated, balanced, bootknife resampling \n');
      else
        fprintf (' Resampling method: Iterated, stratified, balanced, bootknife resampling \n');
      end
    else
      if (isempty (strata))
        fprintf (' Resampling method: Balanced, bootknife resampling \n');
      else
        fprintf (' Resampling method: Stratified, balanced, bootknife resampling \n');
      end
    end
    fprintf (' Number of resamples (outer): %u \n', nboot(1));
    fprintf (' Number of resamples (inner): %u \n', nboot(2));
    if (~ isempty (alpha) && ~ all (isnan (alpha)))
      nalpha = numel (alpha);
      if (nboot(2) > 0)
        if (nalpha > 1)
          fprintf (' Confidence interval (CI) type: Calibrated percentile\n');
        else
          fprintf (' Confidence interval (CI) type: Calibrated percentile (equal-tailed)\n');
        end
      else
        if (nalpha > 1)
          if (strcmpi (bootfun_str, 'mean'))
            fprintf (' Confidence interval (CI) type: Expanded percentile (equal-tailed)\n');
          else
            fprintf (' Confidence interval (CI) type: Percentile (equal-tailed)\n');
          end
        else
          [jnk, warnID] = lastwarn;
          switch warnID
            case 'bootknife:biasfail'
              if (strcmpi (bootfun_str, 'mean'))
                fprintf (' Confidence interval (CI) type: Expanded percentile (equal-tailed)\n');
              else
                fprintf (' Confidence interval (CI) type: Percentile (equal-tailed)\n');
              end
            case 'bootknife:jackfail'
              if (strcmpi (bootfun_str, 'mean'))
                fprintf (' Confidence interval (CI) type: Expanded bias-corrected (BC) \n');
              else
                fprintf (' Confidence interval (CI) type: Bias-corrected (BC) \n');
              end
            otherwise
              if (strcmpi (bootfun_str, 'mean'))
                fprintf (' Confidence interval (CI) type: Expanded bias-corrected and accelerated (BCa) \n');
              else
                fprintf (' Confidence interval (CI) type: Bias-corrected and accelerated (BCa) \n');
              end
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
      if (all (bsxfun (@eq, l, l(1, :))))
        fprintf (' Nominal coverage (and the percentiles used): %.3g%% (%.1f%%, %.1f%%)\n\n', coverage, 100 * l(1, :));
      else
        fprintf (' Nominal coverage: %.3g%%\n\n', coverage);
      end
    end
    fprintf ('Bootstrap Statistics: \n');
    fprintf (' original       bias           std_error      CI_lower       CI_upper    \n');
    for i = 1:m
      fprintf (' %#-+12.6g   %#-+12.6g   %#-+12.6g   %#-+12.6g   %#-+12.6g \n',... 
               [stats.original(i), stats.bias(i), stats.std_error(i), stats.CI_lower(i), stats.CI_upper(i)]);
    end
    fprintf ('\n');
    lastwarn ('', '');  % reset last warning

end

%--------------------------------------------------------------------------

function retval = col2args (func, x, szx)

  % Usage: retval = col2args (func, x, nvar)
  % col2args evaluates func on the columns of x. When nvar > 1, each of the
  % blocks of x are passed to func as a separate arguments. 

  % Extract columns of the matrix into a cell array
  [n, ncols] = size (x);
  xcell = mat2cell (x, n, ncols / sum (szx) * szx);

  % Evaluate column vectors as independent of arguments to bootfun
  retval = func (xcell{:});

end

%--------------------------------------------------------------------------

function [F, x] = empcdf (y)

  % Subfunction to calculate empirical cumulative distribution function

  % Check input argument
  if (~ isa (y, 'numeric'))
    error ('bootknife:empcdf: y must be numeric');
  end
  if (all (size (y) > 1))
    error ('bootknife:empcdf: y must be a vector');
  end
  if (size (y, 2) > 1)
    y = y.';
  end

  % Discard NaN values
  ridx = isnan (y);
  y(ridx) = [];

  % Get size of y
  N = numel (y);

  % Create empirical CDF
  x = sort (y);
  F = linspace (0, 1, N).';

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

  % Apply shrinkage correction
  Y = ((Y - MU) * sqrt (CF)) + MU;

  % Set initial values of X0
  YS = sort (Y, 2);
  X0 = YS(fix ((N - 1) * P) + 1);

  % Perform root finding to get quantiles of the KDE at values of P
  findroot = @(X0, P) fzero (@(X) sum (pnorm (X - Y, 0, BW)) / N - P, X0);
  X = [-Inf, +Inf];
  for i = 1:numel(P)
    if (~ ismember (P(i), [0, 1]))
      X(i) = findroot (X0(i), P(i));
    end
  end

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
  if (exist ('betaincinv', 'file'))
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
      warning ('bootknife:ExpandProbs', ...
          'Could not create studinv function; intervals will not be expanded.');
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
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 95% expanded BCa bootstrap confidence intervals for the mean
%! bootknife (data, 2000, @mean);

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 95% calibrated percentile bootstrap confidence intervals for the mean
%! bootknife (data, [2000, 200], @mean);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 95% calibrated percentile bootstrap confidence intervals for the median
%! ## with smoothing.
%! bootknife (data, [2000, 200], @smoothmedian);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 90% equal-tailed percentile bootstrap confidence intervals for the variance
%! bootknife (data, 2000, {@var, 1}, [0.05, 0.95]);

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 90% BCa bootstrap confidence intervals for the variance
%! bootknife (data, 2000, {@var, 1}, 0.1);

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 90% calibrated equal-tailed percentile bootstrap confidence intervals for
%! ## the variance.
%! bootknife (data, [2000, 200], {@var, 1}, 0.1);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 90% calibrated percentile bootstrap confidence intervals for the variance
%! bootknife (data, [2000, 200], {@var, 1}, [0.05, 0.95]);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input dataset
%! y = randn (20,1); x = randn (20,1); X = [ones(20,1), x];
%!
%! ## 90% BCa confidence interval for regression coefficients 
%! bootknife ({y,X}, 2000, @(y,X) X\y, 0.1); % Could also use @regress


%!demo
%!
%! ## Input bivariate dataset
%! x = [576 635 558 578 666 580 555 661 651 605 653 575 545 572 594].';
%! y = [3.39 3.3 2.81 3.03 3.44 3.07 3 3.43 3.36 3.13 3.12 2.74 2.76 2.88 2.96].'; 
%!
%! ## 95% BCa bootstrap confidence intervals for the correlation coefficient
%! bootknife ({x, y}, 2000, @cor);
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
%! ##      0 33 28 34 4 32 24 47 41 24 26 30 41].';
%! ## boot (1,1,false,1); ci2 = bootknife (A, 20000, {@var,1}, [0.05,0.95]);
%! ## boot (1,1,false,1); ci4 = bootknife (A, 20000, {@var,1}, 0.1);
%! ## boot (1,1,false,1); ci6a = bootknife (A, [20000,500], {@var,1}, 0.1);
%! ## boot (1,1,false,1); ci6b = bootknife (A, [20000,500], {@var,1}, [0.05,0.95]);
%! ##
%! ## Summary of confidence intervals from 'statistics-bootstrap' package for Octave/Matlab
%! ##
%! ## method                               |   0.05 |   0.95 | length | shape |
%! ## -------------------------------------|--------|--------|--------|-------|
%! ## ci2  - percentile (equal-tailed)     |   96.6 |  236.7 |  140.1 |  0.87 |
%! ## ci4  - BCa                           |  115.7 |  266.1 |  150.4 |  1.69 |
%! ## ci6a - calibrated (equal-tailed)     |   82.3 |  256.4 |  174.1 |  0.95 |
%! ## ci6b - calibrated                    |  113.4 |  297.0 |  183.6 |  2.16 |
%! ## -------------------------------------|--------|--------|--------|-------|
%! ## parametric - exact                   |  118.4 |  305.2 |  186.8 |  2.52 |
%! ##
%! ## Simulation results for constructing 90% confidence intervals for the
%! ## variance of a population N(0,1) from 1000 random samples of size 26
%! ## (analagous to the situation above). Simulation performed using the
%! ## bootsim script with nboot of 2000 (for single bootstrap) or [2000,200]
%! ## (for double bootstrap).
%! ##
%! ## method                    | coverage |  lower |  upper | length | shape |
%! ## --------------------------|----------|--------|--------|--------|-------|
%! ## percentile (equal-tailed) |    81.8% |   1.3% |  16.9% |   0.80 |  0.92 |
%! ## BCa                       |    87.3% |   4.2% |   8.5% |   0.86 |  1.85 |
%! ## calibrated (equal-tailed) |    90.5% |   0.5% |   9.0% |   1.06 |  1.06 |
%! ## calibrated                |    90.7% |   5.1% |   4.2% |   1.13 |  2.73 |
%! ## --------------------------|----------|--------|--------|--------|-------|
%! ## parametric - exact        |    90.8% |   3.7% |   5.5% |   0.99 |  2.52 |
%!
%! ## Summary of bias statistics from 'boot' package in R
%! ##
%! ## method             | original |    bias | bias-corrected |
%! ## -------------------|----------|---------|----------------|
%! ## single bootstrap   |   171.53 |   -6.70 |         178.24 |
%! ## double bootstrap   |   171.53 |   -6.83 |         178.36 |
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
%! ## Test for errors when using different functionalities of bootknife
%! ## 'bootknife:parallel' warning turned off in case parallel package is not loaded
%! warning ('off', 'bootknife:parallel')
%! try
%!   y = randn (20,1); 
%!   strata = [1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;3;3;3;3;3];
%!   stats = bootknife (y, 2000, @mean);
%!   stats = bootknife (y, 2000, 'mean');
%!   stats = bootknife (y, 2000, {@var,1});
%!   stats = bootknife (y, 2000, {'var',1});
%!   stats = bootknife (y, 2000, @mean, [], strata);
%!   stats = bootknife (y, 2000, {'var',1}, [], strata);
%!   stats = bootknife (y, 2000, {@var,1}, [], strata, 2);
%!   stats = bootknife (y, 2000, @mean, .1, strata, 2);
%!   stats = bootknife (y, 2000, @mean, [.05,.95], strata, 2);
%!   stats = bootknife (y, [2000,200], @mean, .1, strata, 2);
%!   stats = bootknife (y, [2000,200], @mean, [.05,.95], strata, 2);
%!   stats = bootknife (y(1:5), 2000, @mean, .1);
%!   stats = bootknife (y(1:5), 2000, @mean, [.05,.95]);
%!   stats = bootknife (y(1:5), [2000,200], @mean, .1);
%!   stats = bootknife (y(1:5), [2000,200], @mean, [.05,.95]);
%!   Y = randn (20); 
%!   strata = [1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;3;3;3;3;3];
%!   stats = bootknife (Y, 2000, @mean);
%!   stats = bootknife (Y, 2000, 'mean');
%!   stats = bootknife (Y, 2000, {@var, 1});
%!   stats = bootknife (Y, 2000, {'var',1});
%!   stats = bootknife (Y, 2000, @mean, [], strata);
%!   stats = bootknife (Y, 2000, {'var',1}, [], strata);
%!   stats = bootknife (Y, 2000, {@var,1}, [], strata, 2);
%!   stats = bootknife (Y, 2000, @mean, .1, strata, 2);
%!   stats = bootknife (Y, 2000, @mean, [.05,.95], strata, 2);
%!   stats = bootknife (Y, [2000,200], @mean, .1, strata, 2);
%!   stats = bootknife (Y, [2000,200], @mean, [.05,.95], strata, 2);
%!   stats = bootknife (Y(1:5,:), 2000, @mean, .1);
%!   stats = bootknife (Y(1:5,:), 2000, @mean, [.05,.95]);
%!   stats = bootknife (Y(1:5,:), [2000,200], @mean, .1);
%!   stats = bootknife (Y(1:5,:), [2000,200], @mean, [.05,.95]);
%!   stats = bootknife (Y, 2000, @(Y) mean(Y(:),1)); % Cluster/block resampling
%!   % Y(1,end) = NaN; % Unequal cluster size
%!   %stats = bootknife (Y, 2000, @(Y) mean(Y(:),1,'omitnan'));
%!   y = randn (20,1); x = randn (20,1); X = [ones(20,1), x];
%!   stats = bootknife ({x,y}, 2000, @cor);
%!   stats = bootknife ({x,y}, 2000, @cor, [], strata);
%!   stats = bootknife ({y,x}, 2000, @(y,x) pinv(x)*y); % Could also use @regress
%!   stats = bootknife ({y,X}, 2000, @(y,X) pinv(X)*y);
%!   stats = bootknife ({y,X}, 2000, @(y,X) pinv(X)*y, [], strata);
%!   stats = bootknife ({y,X}, 2000, @(y,X) pinv(X)*y, [], strata, 2);
%!   stats = bootknife ({y,X}, 2000, @(y,X) pinv(X)*y, [.05,.95], strata);
%! catch
%!   warning ('on', 'bootknife:parallel')
%!   rethrow (lasterror)
%! end
%! warning ('on', 'bootknife:parallel')

%!test
%! ## Spatial test data from Table 14.1 of Efron and Tibshirani (1993)
%! ## An Introduction to the Bootstrap in Monographs on Statistics and Applied 
%! ## Probability 57 (Springer)
%! A = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!      0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! ## Nonparametric 90% equal-tailed percentile confidence intervals
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
%! ## Nonparametric 90% calibrated equal-tailed percentile confidence intervals
%! boot (1, 1, false, 1); # Set random seed
%! stats = bootknife(A,[2000,200],{@var,1},0.1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (stats.original, 171.534023668639, 1e-09);
%!   assert (stats.bias, -8.088193809171344, 1e-09);
%!   assert (stats.std_error, 46.53418481731099, 1e-09);
%!   assert (stats.CI_lower, 79.65217027813281, 1e-09);
%!   assert (stats.CI_upper, 260.2974063446341, 1e-09);
%! end
%!
%! ## Nonparametric 90% calibrated percentile confidence intervals
%! boot (1, 1, false, 1); # Set random seed
%! stats = bootknife(A,[2000,200],{@var,1},[0.05,0.95]);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (stats.original, 171.534023668639, 1e-09);
%!   assert (stats.bias, -8.088193809171344, 1e-09);
%!   assert (stats.std_error, 46.53418481731099, 1e-09);
%!   assert (stats.CI_lower, 110.7021156275962, 1e-09);
%!   assert (stats.CI_upper, 305.1908284023669, 1e-09);
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
%! ## Nonparametric 90% equal-tailed percentile confidence intervals
%! ## Percentile intervals on page 266 are 0.524 - 0.928
%! boot (1, 1, false, 1); # Set random seed
%! stats = bootknife({LSAT,GPA},2000,@cor,[0.05,0.95]);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (stats.original, 0.7763744912894071, 1e-09);
%!   assert (stats.bias, -0.008259337758776963, 1e-09);
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
%!   assert (stats.original, 0.7763744912894071, 1e-09);
%!   assert (stats.bias, -0.008259337758776963, 1e-09);
%!   assert (stats.std_error, 0.1420949476115542, 1e-09);
%!   assert (stats.CI_lower, 0.4119228032301614, 1e-09);
%!   assert (stats.CI_upper, 0.9300646701004258, 1e-09);
%! end
%!
%! ## Nonparametric 90% calibrated equal-tailed percentile confidence intervals
%! boot (1, 1, false, 1); # Set random seed
%! stats = bootknife({LSAT,GPA},[2000,500],@cor,0.1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (stats.original, 0.7763744912894071, 1e-09);
%!   assert (stats.bias, -0.00942010836534779, 1e-09);
%!   assert (stats.std_error, 0.1438249935781226, 1e-09);
%!   assert (stats.CI_lower, 0.3730176477191259, 1e-09);
%!   assert (stats.CI_upper, 0.9772521004985222, 1e-09);
%! end
%!
%! ## Nonparametric 90% calibrated percentile confidence intervals
%! boot (1, 1, false, 1); # Set random seed
%! stats = bootknife({LSAT,GPA},[2000,500],@cor,[0.05,0.95]);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (stats.original, 0.7763744912894071, 1e-09);
%!   assert (stats.bias, -0.00942010836534779, 1e-09);
%!   assert (stats.std_error, 0.1438249935781226, 1e-09);
%!   assert (stats.CI_lower, 0.2438194881892977, 1e-09);
%!   assert (stats.CI_upper, 0.944013417640401, 1e-09);
%! end
%! ## Exact intervals based on normal theory are 0.51 - 0.91