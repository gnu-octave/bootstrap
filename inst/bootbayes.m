% -- Function File: bootbayes (y)
% -- Function File: bootbayes (y, X)
% -- Function File: bootbayes (y, X, NBOOT)
% -- Function File: bootbayes (y, X, NBOOT, ALPHA)
% -- Function File: bootbayes (y, X, NBOOT, ALPHA, SEED)
% -- Function File: bootbayes (y, X, NBOOT, ALPHA, SEED, L)
% -- Function File: CI = bootbayes (y, ...)
% -- Function File: [CI, BOOTSTAT] = bootbayes (y, ...)
%
%     'bootbayes (y)' uses Bayesian bootstrap [1] to generate 2000 bootstrap
%     statistics by calculating the mean of the column vector y using weights
%     randomly generated from a symmetric uniform Dirichlet distribution and
%     display the following statistics:
%        • bias: bootstrap estimate of the bias
%        • std_error: bootstrap estimate of the standard error
%        • CI_lower: lower bound of the 95% bootstrap confidence interval
%        • CI_upper: upper bound of the 95% bootstrap confidence interval
%          The confidence intervals, or credible intervals in the context of the
%          Bayesian statistical framework, are equal-tailed percentile intervals.
%
%     'bootbayes (y, X)' specifies the design matrix for least squares
%     regression. X should be a column vector or matrix the same number of
%     rows as y. If the X input argument is empty, the default for X is a
%     is a column of ones (i.e. intercept only) and thus the statistic computed
%     reduces to the mean (as above).
%
%     'bootbayes (y, X, NBOOT)' specifies the number of bootstrap resamples,
%     where NBOOT must be a positive integer. If empty, tHe default value of
%     NBOOT is the scalar: 2000.
%
%     'bootbayes (..., NBOOT, BOOTFUN, ALPHA)' where ALPHA is numeric and
%     sets sets the lower and upper bounds of the confidence interval(s).
%     The value(s) of ALPHA must be between 0 and 1. ALPHA can either be:
%        • scalar: To set the (nominal) central coverage of equal-tailed
%                  percentile confidence intervals to 100*(1-ALPHA)%.
%        • vector: A pair of probabilities defining the (nominal) lower and
%                  upper percentiles of the confidence interval(s) as
%                  100*(ALPHA(1))% and 100*(ALPHA(2))% respectively. 
%        Confidence intervals are not calculated when the value(s) of ALPHA
%        is/are NaN. The default value of  ALPHA is the vector: [.025, .975], 
%        for a 95% confidence interval.
%
%     'bootbayes (..., NBOOT, BOOTFUN, ALPHA, SEED)' initialises the Mersenne
%     Twister random number generator using an integer SEED value so that
%     bootbayes results are reproducible.
%
%     'bootbayes (..., NBOOT, BOOTFUN, ALPHA, SEED, L)' multiplies the
%     regression coefficients by the hypothesis matrix L.
%
%     'CI = bootbayes (STATS, ...) returns the confidence intervals.
%
%     '[CI, BOOTSTAT] = bootbayes (STATS, ...)  also returns BOOTSTAT, a vector
%     or matrix of bootstrap statistics calculated over the bootstrap resamples.
%
%  Bibliography:
%  [1] Rubin (1981) The Bayesian Bootstrap. Ann. Statist. 9(1):130-134
%
%  bootbayes (version 2023.05.02)
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


function [ci, bootstat] = bootbayes (y, X, nboot, alpha, seed, L)

  % Check if running in Octave (else assume Matlab)
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

  % Calculate the length of y
  if (nargin < 1)
    error ('bootbayes: DATA must be provided');
  end
  sz = size (y);
  if ((sz(1) < 2) || (sz (2) > 1))
    error ('bootbayes: y must be a column vector');
  end
  n = numel (y);

  % Evaluate the design matrix
  if (nargin > 1)
    if (isempty (X))
      X = ones (n, 1);
    end
  else
    X = ones (n, 1);
  end

  % Calculate number of parameters
  p = size (X, 2);

  % Evaluate number of bootstrap resamples
  if (nargin > 2)
    if (isempty (nboot))
      nboot = 2000;
    end
    if (~ isa (nboot, 'numeric'))
      error ('bootci: NBOOT must be numeric');
    end
    if (numel (nboot) > 1)
      error ('bootci: NBOOT must be a positive integer');
    end
    if (nboot ~= abs (fix (nboot)))
      error ('bootci: NBOOT must contain positive integers');
    end
  else
    nboot = 2000;
  end

  % Evaluate alpha
  if ((nargin < 4) || isempty (alpha))
    alpha = [0.025, 0.975];
    nalpha = 2;
  else
    nalpha = numel (alpha);
    if (~ isa (alpha, 'numeric') || (nalpha > 2))
      error ('bootbayes: ALPHA must be a scalar (two-tailed probability) or a vector (pair of probabilities)');
    end
    if (size (alpha, 1) > 1)
      alpha = alpha.';
    end
    if (any ((alpha < 0) | (alpha > 1)))
      error ('bootbayes: Value(s) in ALPHA must be between 0 and 1');
    end
    if (nalpha > 1)
      % alpha is a pair of probabilities
      % Make sure probabilities are in the correct order
      if (alpha(1) > alpha(2) )
        error ('bootbayes: The pair of probabilities must be in ascending numeric order');
      end
    end
  end
  if (nalpha < 2)
    % Create equal-tailed probabilities for the percentiles
    l = [alpha / 2, 1 - alpha / 2];
  else
    l = alpha;
  end

  % Set random seed
  if (nargin > 4)
    if (~ isempty (seed))
      if (ISOCTAVE)
        rande ('seed', seed);
      else
        rng ('default');
      end
    end
  end

  % Evaluate hypothesis matrix (L)
  if (nargin < 6)
    % If L is not provided, set L to unity
    L = 1;
  else
    if (isempty (H))
      L = 1;
    end
    % Calculate number of parameters
    p = size (L, 1);
  end

  % Create weighted least squares anonymous function
  bootfun = @(w) lmfit (X, y, diag (w), L);

  % Calculate estimate(s)
  T0 = bootfun (ones (n, 1));

  % Create weights by randomly sample from a symmetric uniform Dirichlet distribution
  r = exprnd (1, n, nboot); % Equal to gamma distribution with scale = shape = 1
  W = bsxfun (@rdivide, r, sum (r));

  % Compute bootstat
  bootstat = cell2mat (cellfun (bootfun, num2cell (W, 1), 'UniformOutput', false));

  % Bootstrap bias estimation
  bias = mean (bootstat, 2) - T0;

  % Bootstrap standard error
  se = std (bootstat, 0, 2);

  % Compute confidence intervals
  ci = nan (p, 2);
  if (~ isnan (alpha))
    ci = zeros (p, 2);
    for j = 1:p
      [cdf, t1] = empcdf (bootstat(j, :));
      ci(j, :) = arrayfun (@(p) interp1 (cdf, t1, p, 'linear'), l);
    end
  end

  % Prepare output arguments
  stats = struct;
  stats.original = T0;
  stats.bias = bias;
  stats.std_error = se;
  stats.CI_lower = ci(:, 1);
  stats.CI_upper = ci(:, 2);

  % Print output if no output arguments are requested
  if (nargout == 0) 
    print_output (stats, nboot, alpha, l, p);
  end

end

%--------------------------------------------------------------------------

%% FUNCTION TO FIT THE LINEAR MODEL

function b = lmfit (X, y, W, L)

  % Get model coefficients by solving the linear equation by matrix arithmetic
  % If optional arument W is provided, it should be a diagonal matrix of
  % weights or a positive definite covariance matrix
  if (nargin < 3)
    % If no weights are provided, create an identity matrix
    n = numel (y);
    W = eye (n);
  end
  if (nargin < 4)
    % If no hypothesis matrix is provided, set H to unity
    L = 1;
  end
  
  % Solve linear equation to minimize weighted least squares
  b = L * pinv (X' * W * X) * (X' * W * y);

end

%--------------------------------------------------------------------------

%% FUNCTION TO OBTAIN EMPIRICAL CUMULATIVE DISTRIBUTION FUNCTION

function [F, x] = empcdf (y)

  % Subfunction to calculate empirical cumulative distribution function

  % Check input argument
  if (~ isa (y, 'numeric'))
    error ('bootbayes:empcdf: y must be numeric');
  end
  if (all (size (y) > 1))
    error ('bootbayes:empcdf: y must be a vector');
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

%% FUNCTION TO PRINT OUTPUT

function print_output (stats, nboot, alpha, l, p)

    fprintf (['\nSummary of Bayesian bootstrap estimates of bias and precision for linear models\n',...
              '*******************************************************************************\n\n']);
    fprintf ('Bootstrap settings: \n');
    fprintf (' Function: (X'' * W * y) / (X'' * W * X)\n');
    fprintf (' Resampling method: Bayesian bootstrap (symmetric uniform Dirichlet)\n')
    fprintf (' Number of resamples: %u \n', nboot)
    if (~ isempty (alpha) && ~ all (isnan (alpha)))
      nalpha = numel (alpha);
      fprintf (' Confidence interval (CI) type: Percentile (equal-tailed)\n');
      if (nalpha > 1)
        % alpha is a vector of probabilities
        coverage = 100 * abs (alpha(2) - alpha(1));
      else
        % alpha is a two-tailed probability
        coverage = 100 * (1 - alpha);
      end
      fprintf (' Nominal coverage (and the percentiles used): %.3g%% (%.1f%%, %.1f%%)\n\n', coverage, 100 * l);
    end
    fprintf ('Bootstrap Statistics: \n');
    fprintf (' original       bias           std_error      CI_lower       CI_upper    \n');
    for j = 1:p
      fprintf (' %#-+12.6g   %#-+12.6g   %#-+12.6g   %#-+12.6g   %#-+12.6g \n',... 
               [stats.original(j), stats.bias(j), stats.std_error(j), stats.CI_lower(j), stats.CI_upper(j)]);
    end
    fprintf ('\n');

end

%--------------------------------------------------------------------------

%!demo
%!
%! ## Input bivariate dataset
%! heights = [183, 192, 182, 183, 177, 185, 188, 188, 182, 185].';
%!
%! ## 95% bootstrap confidence interval for the mean 
%! bootbayes(heights);
%!
%! ## Please be patient, the calculations will be completed soon...
