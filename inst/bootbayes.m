% -- Function File: bootbayes (y)
% -- Function File: bootbayes (y, X)
% -- Function File: bootbayes (y, X, NBOOT)
% -- Function File: bootbayes (y, X, NBOOT, PROB)
% -- Function File: bootbayes (y, X, NBOOT, PROB, PRIOR)
% -- Function File: bootbayes (y, X, NBOOT, PROB, PRIOR, SEED)
% -- Function File: bootbayes (y, X, NBOOT, PROB, PRIOR, SEED, L)
% -- Function File: STATS = bootbayes (y, ...)
% -- Function File: [STATS, BOOTSTAT] = bootbayes (y, ...)
%
%     'bootbayes (y)' performs Bayesian nonparametric bootstrap [1] to create
%     2000 bootstrap statistics, each representing the weighted mean of the
%     column vector, y, using a vector of weights randomly generated from a
%     symmetric Dirichlet distribution. The resulting bootstrap (or posterior
%     [1,2]) distribution(s) is/are summarised with the following statistics
%     printed to the standard output:
%        • original: the mean of the data vector y
%        • bias: bootstrap estimate(s) of the bias
%        • median: the median of the posterior distribution(s)
%        • CI_lower: lower bound(s) of the 95% credible interval
%        • CI_upper: upper bound(s) of the 95% credible interval
%          By default, the credible intervals are shortest probability intervals,
%          which represent a more computationally stable version of the highest
%          posterior density interval [3].
%
%     'bootbayes (y, X)' also specifies the design matrix (X) for least squares
%     regression of y on X. X should be a column vector or matrix the same
%     number of rows as y. If the X input argument is empty, the default for X
%     is a column of ones (i.e. intercept only) and thus the statistic computed
%     reduces to the mean (as above). The statistics calculated and returned in
%     the output then relate to the coefficients from the regression of y on X.
%
%     'bootbayes (y, X, NBOOT)' specifies the number of bootstrap resamples,
%     where NBOOT must be a positive integer. If empty, the default value of
%     NBOOT is 2000.
%
%     'bootbayes (..., NBOOT, PROB)' where PROB is numeric and sets the lower
%     lower and upper bounds of the credible interval(s). The value(s) of
%     PROB must be between 0 and 1. PROB can either be:
%        • scalar: To set the central mass of shortest probability intervals
%                  (SPI) to 100*(1-PROB)%
%        • vector: A pair of probabilities defining the lower and upper
%                  percentiles of the credible interval(s) as 100*(PROB(1))%
%                  and 100*(PROB(2))% respectively. 
%          Credible intervals are not calculated when the value(s) of PROB
%          is/are NaN. The default value of PROB is 0.95.
%
%     'bootbayes (..., NBOOT, PROB, PRIOR)' accepts a positive real numeric
%     scalar to parametrize the form of the symmetric Dirichlet distribution.
%     The Dirichlet distribution is the conjugate PRIOR used to randomly
%     generate weights for linear least squares fitting to the observed data,
%     and subsequently to estimate the posterior for the regression coefficients
%     by Bayesian bootstrap. If PRIOR is not provided, or is empty, it will be
%     set to 1, corresponding to Bayes rule: a uniform (or flat) Dirichlet
%     distribution (in the range [0, 1]). For a stronger prior, set PRIOR to 
%     > 1, for example, if the experiment is confirmatory and reproduced. For a
%     weaker prior, set PRIOR to < 1 (e.g. 0.5 for Jeffrey's prior). Jeffrey's
%     prior may be appropriate for estimates from small samples (n < 10) in 
%     exploratory settings, where the amount of data may be assumed to
%     inadequately define the parameter space. 
%
%     'bootbayes (..., NBOOT, PROB, PRIOR, SEED)' initialises the Mersenne
%     Twister random number generator using an integer SEED value so that
%     'bootbayes' results are reproducible.
%
%     'bootbayes (..., NBOOT, PROB, PRIOR, SEED, L)' multiplies the regression
%     coefficients by the hypothesis matrix L. If L is not provided or is empty,
%     it will assume the default value of 1.
%
%     'STATS = bootbayes (STATS, ...) returns a structure with the following
%     fields (defined above): original, bias, median, CI_lower, CI_upper.
%
%     '[STATS, BOOTSTAT] = bootbayes (STATS, ...)  also returns the a vector (or
%     matrix) of bootstrap statistics (BOOTSTAT) calculated over the bootstrap
%     resamples.
%
%  Bibliography:
%  [1] Rubin (1981) The Bayesian Bootstrap. Ann. Statist. 9(1):130-134
%  [2] Weng (1989) On a Second-Order Asymptotic property of the Bayesian
%        Bootstrap Mean. Ann. Statist. 17(2):705-710
%  [3] Liu, Gelman & Zheng (2015). Simulation-efficient shortest probability
%        intervals. Statistics and Computing, 25(4), 809–819. 
%
%  bootbayes (version 2023.05.14)
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


function [stats, bootstat] = bootbayes (y, X, nboot, prob, prior, seed, L)

  % Check the number of function arguments
  if (nargin < 1)
    error ('bootknife: y must be provided');
  end
  if (nargin > 7)
    error ('bootbayes: Too many input arguments')
  end
  if (nargout > 2)
    error ('bootbayes: Too many output arguments')
  end

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
      error ('bootbayes: NBOOT must be numeric');
    end
    if (numel (nboot) > 1)
      error ('bootbayes: NBOOT must be a positive integer');
    end
    if (nboot ~= abs (fix (nboot)))
      error ('bootbayes: NBOOT must contain positive integers');
    end
  else
    nboot = 2000;
  end

  % Evaluate prob
  if ((nargin < 4) || isempty (prob))
    prob = 0.95;
    nprob = 1;
  else
    nprob = numel (prob);
    if (~ isa (prob, 'numeric') || (nprob > 2))
      error ('bootbayes: PROB must be a scalar or a vector of length 2');
    end
    if (size (prob, 1) > 1)
      prob = prob.';
    end
    if (any ((prob < 0) | (prob > 1)))
      error ('bootbayes: Value(s) in PROB must be between 0 and 1');
    end
    if (nprob > 1)
      % PROB is a pair of probabilities
      % Make sure probabilities are in the correct order
      if (prob(1) > prob(2) )
        error ('bootbayes: The pair of probabilities must be in ascending numeric order');
      end
    end
  end

  % Evaluate or set prior
  % Set the prior based on our prior expectation that, depending on the sample
  % size, sample variance underestimates the population variance 
  if (nargin < 5)
    prior = 1; % Bayes flat/uniform prior
  else
    if (isempty (prior))
      prior = 1; % Bayes flat/uniform prior
    end
    if (~ isa (prior, 'numeric'))
      error ('bootbayes: PRIOR must be numeric');
    end
    if (numel (prior) > 1)
      error ('bootbayes: PRIOR must be scalar');
    end
    if (prior ~= abs (prior))
      error ('bootbayes: PRIOR must be positive');
    end
    if ~ (prior > 0)
      error ('bootbayes: PRIOR must be greater than zero')
    end
  end

  % Set random seed
  if (nargin > 5)
    if (~ isempty (seed))
      if (ISOCTAVE)
        randg ('seed', seed);
      else
        rng ('default');
      end
    end
  end

  % Evaluate hypothesis matrix (L)
  if (nargin < 7)
    % If L is not provided, set L to 1
    L = 1;
  else
    if (isempty (L))
      L = 1;
    end
    % Calculate number of parameters
    p = size (L, 1);
  end

  % Create weighted least squares anonymous function
  bootfun = @(w) lmfit (X, y, diag (w), L);

  % Calculate estimate(s)
  original = bootfun (ones (n, 1));

  % Create weights by randomly sampling from a symmetric Dirichlet distribution.
  % This can be achieved by normalizing a set of randomly generated values from
  % a Gamma distribution to their sum.
  if (ISOCTAVE)
    r = randg (prior, n, nboot);
  else
    r = gamrnd (prior, 1, n, nboot);
  end
  W = bsxfun (@rdivide, r, sum (r));

  % Compute bootstat
  bootstat = cell2mat (cellfun (bootfun, num2cell (W, 1), 'UniformOutput', false));

  % Bootstrap bias estimation
  bias = mean (bootstat, 2) - original;

  % Compute credible intervals
  % https://discourse.mc-stan.org/t/shortest-posterior-intervals/16281/16
  ci = nan (p, 2);
  bootstat = sort (bootstat, 2);
  gap = round (prob * nboot);
  for j = 1:p
    if (nprob > 1)
      % Percentile intervals
      if (~ isnan (prob))
        ci(j, :) = bootstat(j, gap);
      end
    else
      % Shortest probability interval
      width = bootstat(j, (gap + 1) : nboot) - bootstat(j, 1 : (nboot - gap));
      index = min (find (width == min (width)));
      if (~ isnan (prob))
        ci(j, :) = bootstat(j, [index, index + gap]);
      end
    end
  end
  
  % Prepare output arguments
  stats = struct;
  stats.original = original;
  stats.bias = bias;
  stats.median = median (bootstat, 2);
  stats.CI_lower = ci(:, 1);
  stats.CI_upper = ci(:, 2);

  % Print output if no output arguments are requested
  if (nargout == 0) 
    print_output (stats, nboot, prob, prior, p, L);
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
    % If no hypothesis matrix (L) is provided, set L to 1
    L = 1;
  end
  
  % Solve linear equation to minimize weighted least squares
  b = L * pinv (X' * W * X) * (X' * W * y);

end

%--------------------------------------------------------------------------

%% FUNCTION TO PRINT OUTPUT

function print_output (stats, nboot, prob, prior, p, L)

    fprintf (['\nSummary of Bayesian bootstrap estimates of bias and precision for linear models\n',...
              '*******************************************************************************\n\n']);
    fprintf ('Bootstrap settings: \n');
    fprintf (' Function: L * pinv (X'' * W * X) * (X'' * W * y)\n');
    fprintf (' Resampling method: Bayesian bootstrap\n')
    fprintf (' Prior: Symmetric Dirichlet distribution (a = %.3g)\n', prior)
    fprintf (' Number of resamples: %u \n', nboot)
    if (~ isempty (prob) && ~ all (isnan (prob)))
      nprob = numel (prob);
      if (nprob > 1)
        % prob is a vector of probabilities
        fprintf (' Credible interval (CI) type: Percentile interval\n');
        mass = 100 * abs (prob(2) - prob(1));
        fprintf (' Credible interval: %.3g%% (%.1f%%, %.1f%%)\n', mass, 100 * prob);
      else
        % prob is a two-tailed probability
        fprintf (' Credible interval (CI) type: Shortest probability interval\n');
        mass = 100 * prob;
        fprintf (' Credible interval: %.3g%%\n', mass);
      end
    end
    fprintf ('\nPosterior Statistics: \n');
    fprintf (' original       bias           median         CI_lower       CI_upper\n');
    for j = 1:p
      fprintf (' %#-+12.6g   %#-+12.6g   %#-+12.6g   %#-+12.6g   %#-+12.6g \n',... 
                 [stats.original(j), stats.bias(j), stats.median(j), stats.CI_lower(j), stats.CI_upper(j)]);
    end
    fprintf ('\n');

end

%--------------------------------------------------------------------------

%!demo
%!
%! ## Input univariate dataset
%! heights = [183, 192, 182, 183, 177, 185, 188, 188, 182, 185].';
%!
%! ## 95% credible interval for the mean 
%! bootbayes(heights);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input bivariate dataset
%! X = [ones(43,1),...
%!     [01,02,03,04,05,06,07,08,09,10,11,...
%!      12,13,14,15,16,17,18,19,20,21,22,...
%!      23,25,26,27,28,29,30,31,32,33,34,...
%!      35,36,37,38,39,40,41,42,43,44]'];
%! y = [188.0,170.0,189.0,163.0,183.0,171.0,185.0,168.0,173.0,183.0,173.0,...
%!     173.0,175.0,178.0,183.0,192.4,178.0,173.0,174.0,183.0,188.0,180.0,...
%!     168.0,170.0,178.0,182.0,180.0,183.0,178.0,182.0,188.0,175.0,179.0,...
%!     183.0,192.0,182.0,183.0,177.0,185.0,188.0,188.0,182.0,185.0]';
%!
%! ## 95% credible interval for the regression coefficents
%! bootbayes(y,X,1000);
%!
%! ## Please be patient, the calculations will be completed soon...

%!test
%! ## Test calculations of statistics for the mean
%!
%! ## Input univariate dataset
%! heights = [183, 192, 182, 183, 177, 185, 188, 188, 182, 185].';
%!
%! ## 95% credible interval for the mean 
%! stats = bootbayes(heights);
%! stats = bootbayes(heights,[]);
%! stats = bootbayes(heights,[],2000);
%! stats = bootbayes(heights,[],2000,0.05);
%! stats = bootbayes(heights,[],2000,[0.025,0.975]);
%! stats = bootbayes(heights,[],[],[]);
%! [stats,bootstat] = bootbayes(heights);

%!test
%! ## Test calculations of statistics for linear regression
%!
%! ## Input bivariate dataset
%! X = [ones(43,1),...
%!     [01,02,03,04,05,06,07,08,09,10,11,...
%!      12,13,14,15,16,17,18,19,20,21,22,...
%!      23,25,26,27,28,29,30,31,32,33,34,...
%!      35,36,37,38,39,40,41,42,43,44]'];
%! y = [188.0,170.0,189.0,163.0,183.0,171.0,185.0,168.0,173.0,183.0,173.0,...
%!     173.0,175.0,178.0,183.0,192.4,178.0,173.0,174.0,183.0,188.0,180.0,...
%!     168.0,170.0,178.0,182.0,180.0,183.0,178.0,182.0,188.0,175.0,179.0,...
%!     183.0,192.0,182.0,183.0,177.0,185.0,188.0,188.0,182.0,185.0]';
%!
%! ## 95% credible interval for the mean 
%! stats = bootbayes(y,X);
%! stats = bootbayes(y,X,2000);
%! stats = bootbayes(y,X,2000,0.05);
%! stats = bootbayes(y,X,2000,[0.025,0.975]);
%! stats = bootbayes(y,X,[],[]);
%! [stats,bootstat] = bootbayes(y,X);