% -- Function File: bootclust (DATA, CLUSTID)
% -- Function File: bootclust ({D1, D2, ...}, CLUSTID)
% -- Function File: bootclust (..., CLUSTID, BOOTFUN)
% -- Function File: bootclust (..., CLUSTID, {BOOTFUN, ...})
% -- Function File: bootclust (..., CLUSTID, ..., NBOOT)
% -- Function File: bootclust (..., CLUSTID, ..., NBOOT, ALPHA)
% -- Function File: bootclust (..., CLUSTID, ..., NBOOT, ALPHA, SEED)
% -- Function File: STATS = bootclust (...)
% -- Function File: [STATS, BOOTSTAT] = bootclust (...)
%
%     'bootclust (DATA, CLUSTID)' uses a variant of nonparametric bootstrap,
%     called bootknife [1], to generate 1999 resamples from clusters of the
%     rows of DATA (column vector or matrix) defined by CLUSTID and compute
%     their means and display the following statistics:
%        - original: the original estimate(s) calculated by BOOTFUN and the DATA
%        - bias: bootstrap bias of the estimate(s)
%        - std_error: bootstrap standard error of the estimate(s)
%        - CI_lower: lower bound(s) of the 95% bootstrap confidence interval
%        - CI_upper: upper bound(s) of the 95% bootstrap confidence interval
%        CLUSTID must be a vector or cell array of numbers or strings
%        respectively to be used as cluster labels or identifiers. Rows in
%        DATA with the same CLUSTID value are treated as clusters of
%        observations that are resampled together. Since this function uses
%        bootknife resampling [1], the variance of the bootstrap statistics is
%        an unbiased estimator of the sampling variance, and so the standard
%        errors are also unbiased. The confidence intervals returned are
%        percentile intervals computed from a kernel density estimate of the
%        bootstrap statistics (with shrinkage correction). These intervals
%        resemble the normal intervals (i.e. original +/- std_error * norminv
%        (1-alpha/2)) at very small sample sizes, but become increasingly more
%        like percentile intervals as a function of increasing sample size.
%
%     'bootclust ({D1, D2, ...}, CLUSTID)' resamples from the clusters of the
%     data vectors D1, D2 etc and the resamples are passed onto BOOTFUN as
%     multiple data input arguments. All data vectors and matrices (D1, D2 etc)
%     must have the same number of rows.
%
%     'bootclust (..., CLUSTID, BOOTFUN)' also specifies BOOTFUN: the function
%     calculated on the original sample and the bootstrap resamples. BOOTFUN
%     must be either a:
%       <> function handle,
%       <> string of function name, or
%       <> a cell array where the first cell is one of the above function
%          definitions and the remaining cells are (additional) input arguments 
%          to that function (other than the data arguments).
%        In all cases BOOTFUN must take DATA for the initial input argument(s).
%        BOOTFUN can return a scalar or any multidimensional numeric variable,
%        but the output will be reshaped as a column vector. BOOTFUN must
%        calculate a statistic representative of the finite data sample; it
%        should NOT be an estimate of a population parameter (unless they are
%        one of the same). If BOOTFUN is @mean or 'mean', narrowness bias of
%        the confidence intervals for single bootstrap are reduced by expanding
%        the probabilities of the percentiles using Student's t-distribution
%        [2]. By default, BOOTFUN is @mean.
%
%     'bootclust (..., CLUSTID, BOOTFUN, NBOOT)' specifies the number of
%     bootstrap resamples, where NBOOT is a scalar, positive integer
%     corresponding to the number of bootstrap resamples. THe default value
%     of NBOOT is the scalar: 1999.
%
%     'bootclust (..., CLUSTID, BOOTFUN, NBOOT, ALPHA)', where ALPHA is numeric
%     and sets the lower and upper bounds of the confidence interval(s). The
%     value(s) of ALPHA must be between 0 and 1. ALPHA can either be:
%       <> scalar: To set the (nominal) central coverage of equal-tailed
%                  percentile confidence intervals to 100*(1-ALPHA)%.
%       <> vector: A pair of probabilities defining the (nominal) lower and
%                  upper percentiles of the confidence interval(s) as
%                  100*(ALPHA(1))% and 100*(ALPHA(2))% respectively.
%        The default value of ALPHA is the vector: 0.05 for an equitailed 95%
%        percentile confidence interval.
%
%     'bootclust (..., CLUSTID, BOOTFUN, NBOOT, ALPHA, SEED)' initialises the
%     Mersenne Twister random number generator using an integer SEED value so
%     that bootci results are reproducible.
%
%     'STATS = bootclust (...)' returns a structure with the following fields
%     (defined above): original, bias, std_error, CI_lower, CI_upper.
%
%     '[STATS, BOOTSTAT] = bootclust (...)' returns BOOTSTAT, a vector or matrix
%     of bootstrap statistics calculated over the bootstrap resamples.
%
%  BIBLIOGRAPHY:
%  [1] Hesterberg T.C. (2004) Unbiasing the Bootstrap—Bootknife Sampling 
%        vs. Smoothing; Proceedings of the Section on Statistics & the 
%        Environment. Alexandria, VA: American Statistical Association.
%  [2] Hesterberg, Tim (2014), What Teachers Should Know about the 
%        Bootstrap: Resampling in the Undergraduate Statistics Curriculum, 
%        http://arxiv.org/abs/1411.5279
%
%  bootclust (version 2023.09.18)
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

function [stats, bootstat] = bootclust (x, clustid, bootfun, nboot, alpha, seed)

  % Check if we are running Octave or Matlab
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));
  % Store local functions in a stucture for parallel processes
  localfunc = struct ('col2args', @col2args, ...
                      'kdeinv', @kdeinv, ...
                      'ExpandProbs', @ExpandProbs);

  % Check if we have parallel processing capabilities
  PARALLEL = false; % Default
  if (ISOCTAVE)
    software = pkg ('list');
    names = cellfun (@(S) S.name, software, 'UniformOutput', false);
    status = cellfun (@(S) S.loaded, software, 'UniformOutput', false);
    index = find (~ cellfun (@isempty, regexpi (names, '^parallel')));
    if ( (~ isempty (index)) && (logical (status{index})) )
      PARALLEL = true;
    end
  else
    try 
      pool = gcp ('nocreate'); 
      PARALLEL = ~ isempty (pool);
    catch
      % Do nothing
    end
  end

  % Check the number of function arguments
  if (nargin < 1)
    error ('bootclust: DATA must be provided');
  end
  if (nargin > 6)
    error ('bootclust: Too many input arguments')
  end
  if (nargout > 2)
    error ('bootclust: Too many output arguments')
  end

  % BOOTFUN input argument
  if ((nargin < 2) || isempty (bootfun))
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
      error ('bootclust: BOOTFUN must be a function name or function handle')
    end
  end

  % NBOOT input argument
  if ((nargin < 3) || isempty (nboot))
    nboot = 1999;
  else
    if (~ isa (nboot, 'numeric'))
      error ('bootclust: NBOOT must be numeric');
    end
    if (numel (nboot) > 1)
      error ('bootclust: NBOOT cannot contain more than 1 value');
    end
    if (nboot ~= abs (fix (nboot)))
      error ('bootclust: NBOOT must contain positive integers');
    end    
  end
  if (~ all (size (nboot) == [1, 1]))
    error ('bootclust: NBOOT must be a scalar value')
  end

  % ALPHA input argument
  if ( (nargin < 5) || isempty (alpha) )
    alpha = .05;
  end
  nalpha = numel (alpha);
  if (~ isa (alpha, 'numeric') || (nalpha > 2))
    error (cat (2, 'bootclust: ALPHA must be a scalar (two-tailed', ...
                   'probability) or a vector (pair of probabilities)'))
  end
  if (size (alpha, 1) > 1)
    alpha = alpha.';
  end
  if (any (isnan (alpha)))
    error ('bootclust: ALPHA cannot contain NaN values');
  end
  if (any ((alpha < 0) | (alpha > 1)))
    error ('bootclust: Value(s) in ALPHA must be between 0 and 1');
  end
  if (nalpha > 1)
    % alpha is a pair of probabilities
    % Make sure probabilities are in the correct order
    if (alpha(1) > alpha(2) )
      error (cat (2, 'bootknife: The pair of probabilities must be', ...
                     ' in ascending numeric order'))
    end
    probs = alpha;
    alpha = 1 - probs(2) + probs(1);
  else
    probs = [alpha / 2 , 1 - alpha / 2];
  end

  % Initialise the random number generator with the SEED (if provided)
  if ( (nargin > 5) && (~ isempty (seed)) )
    boot (1, 1, true, seed);
  end

  % If DATA is a cell array of equal size colunmn vectors, convert the cell
  % array to a matrix and redefine bootfun to parse multiple input arguments
  if (iscell (x))
    szx = cellfun (@(x) size (x, 2), x);
    x = [x{:}];
    bootfun = @(x) localfunc.col2args (bootfun, x, szx);
  else
    szx = size (x, 2);
  end

  % Determine properties of the DATA (x)
  [n, nvar] = size (x);
  if (n < 2)
    error ('bootclust: DATA must be numeric and contain > 1 row')
  end

  % Sort rows of CLUSTID and the DATA
  if ((nargin < 2) || isempty (clustid))
    clustid = (1:n)';
  else
    if ( any (size (clustid) ~= [n, 1]) )
      error (cat (2, 'bootclust: CLUSTID must be a column vector with', ...
                     ' the same number of rows as DATA'))
    end
    [clustid, idx] = sort (clustid);
    x = x(idx,:);
  end

  % Evaluate definition of the sampling units (e.g. clusters) of x 
  [ux, jnk, ic] = unique (clustid);
  nx = numel (ux);

  % Calculate the number of elements in the return value of bootfun and check
  % whether function evaluations can be vectorized
  T0 = bootfun (x);
  m = numel (T0);
  if (nvar > 1)
    M = cell2mat (cellfun (@(i) repmat (x(:, i), 1, 2), ... 
                 num2cell (1 : nvar), 'UniformOutput', false));
  else 
    M = repmat (x, 1, 2); 
  end
  if (any (szx > 1))
    VECTORIZED = false;
  else
    try
      chk = bootfun (M);
      if (all (size (chk) == [size(T0, 1), 2]) && all (chk == bootfun (x)))
        VECTORIZED = true;
      else
        VECTORIZED = false;
      end
    catch
      VECTORIZED = false;
    end
  end
  if (m > 1)
    % Vectorized along the dimension of the return values of bootfun so
    % reshape the output to be a column vector before proceeding with bootstrap
    if (size (T0, 2) > 1)
      bootfun = @(x) reshape (bootfun (x), [], 1);
      T0 = reshape (T0, [], 1);
      VECTORIZED = false;
    end
  end
  % Check if we can vectorize function evaluations
  if (any (diff (accumarray (ic, 1))))
    VECTORIZED = false;
  end

  % Convert x to a cell array of clusters
  x = mat2cell (x, accumarray (ic, 1));

  % Perform resampling of clusters
  bootsam = boot (nx, nboot, true);
  X = arrayfun (@(i) x(bootsam(i, :)), 1 : nx, 'UniformOutput', false);
  X = [X{:}]';

  % Perform the function evaluations
  if VECTORIZED
    X = reshape (vertcat (X{:}), n, nboot * nvar);
    bootstat = bootfun (X);
  else
    if (PARALLEL)
      if (ISOCTAVE)
        % Set unique random seed for each parallel thread
        bootstat = cell2mat (pararrayfun (nproc, ...
                                          @(b) bootfun (vertcat (X{:,b})), ...
                                          1 : nboot, 'UniformOutput', false));
      else
        % Set unique random seed for each parallel thread
        bootstat = zeros (m, nboot);
        parfor b = 1:nboot
          bootstat(:, b) = bootfun (vertcat (X{:,b}));
        end
      end
    else
      bootstat = cell2mat (arrayfun (@(b) bootfun (vertcat (X{:,b})), ...
                                          1 : nboot, 'UniformOutput', false));
    end
  end

  % Remove bootstrap statistics that contain NaN or inf
  ridx = any (or (isnan (bootstat), isinf (bootstat)) , 1);
  bootstat(:, ridx) = [];
  if (isempty (bootstat))
    error ('bootclust: BOOTFUN returned NaN or inf for all bootstrap resamples')
  end
  nboot = nboot - sum (ridx);

  % Bootstrap bias estimation
  bias = mean (bootstat, 2) - T0;

  % Bootstrap standard error
  se = std (bootstat, 0, 2);  % Unbiased since we used bootknife resampling

  % If bootfun is the arithmetic meam, expand the probabilities of the 
  % percentiles for the confidence intervals using Student's t-distribution
  if (strcmpi (bootfun_str, 'mean'))
    probs = localfunc.ExpandProbs (probs, nx);
  end

  % Intervals constructed from kernel density estimate of the bootstrap
  % statistics (with shrinkage correction)
  ci = nan (m, 2);
  for j = 1:m
    try
      ci(j, :) = localfunc.kdeinv (probs, bootstat(j, :), ...
                               se(j) * sqrt (1 / (nx - 1)), 1 - 1 / (nx - 1));
    catch
      % Linear interpolation (legacy)
      fprintf (strcat ('Note: Falling back to linear interpolation to', ...
                           ' calculate percentiles for interval pair %u\n'), j);
      [t1, cdf] = bootcdf (bootstat(j, :), true, 1);
      ci(j, 1) = interp1 (cdf, t1, probs(1), 'linear', min (t1));
      ci(j, 2) = interp1 (cdf, t1, probs(2), 'linear', max (t1));
    end
  end

  % Create STATS output structure
  stats = struct;
  stats.original = T0;
  stats.bias = bias;          % Bootstrap bias estimation
  stats.std_error = se;       % Bootstrap standard error
  stats.CI_lower = ci(:, 1);  % Lower percentile
  stats.CI_upper = ci(:, 2);  % Upper percentile

  % Print output if no output arguments are requested
  if (nargout == 0) 
    print_output (stats, nboot, alpha, probs, m, bootfun_str);
  end

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

  % Modify ALPHA to adjust tail probabilities assuming that the kurtosis
  % of the sampling distribution scales with degrees of freedom like the
  % t-distribution. This is related in concept to ExpandProbs in the
  % R package 'resample':
  % www.rdocumentation.org/packages/resample/versions/0.6/topics/ExpandProbs

  % Get size of P
  sz = size (P);

  % Create required distribution functions
  stdnormcdf = @(X) 0.5 * (1 + erf (X / sqrt (2)));
  stdnorminv = @(P) sqrt (2) * erfinv (2 * P - 1);
  if (exist ('betaincinv', 'file'))
    studinv = @(P, DF) sign (P - 0.5) * ...
                sqrt ( DF ./ betaincinv (2 * min (P, 1 - P), DF / 2, 0.5) - DF);
  else
    % Earlier versions of Matlab do not have betaincinv
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

function print_output (stats, nboot, alpha, probs, m, bootfun_str)

    fprintf (cat (2, '\nSummary of cluster bootstrap estimates of', ...
                     ' bias and precision\n', ...
                     '*************************************************', ...
                     '*****************************\n\n'));
    fprintf ('Bootstrap settings: \n');
    fprintf (' Function: %s\n', bootfun_str);
    fprintf (' Resampling method: Balanced, bootknife resampling \n');
    fprintf (' Number of resamples: %u \n', nboot(1));
    if (strcmpi (bootfun_str, 'mean'))
      fprintf (' Confidence interval (CI) type: Expanded percentile\n');
    else
      fprintf (' Confidence interval (CI) type: Percentile\n');
    end
    coverage = 100 * (1 - alpha);
    fprintf (cat (2, ' Nominal coverage (and the percentiles used):', ...
                      ' %.3g%% (%.1f%%, %.1f%%)\n\n'), coverage, 100 * probs);
    fprintf ('Bootstrap Statistics: \n');
    fprintf (cat (2, ' original     bias         std_error    CI_lower', ...
                     '     CI_upper  \n'));
    for i = 1:m
      fprintf (cat (2, ' %#-+10.4g   %#-+10.4g   %#-+10.4g   %#-+10.4g', ...
                     '   %#-+10.4g \n'), [stats.original(i), stats.bias(i), ...
                     stats.std_error(i), stats.CI_lower(i), stats.CI_upper(i)]);
    end
    fprintf ('\n');

end

%--------------------------------------------------------------------------