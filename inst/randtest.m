% Performs permutation or randomization tests for regression coefficients.
%
% -- Function File: PVAL = randtest (X, Y)
% -- Function File: PVAL = randtest (X, Y, NREPS)
% -- Function File: PVAL = randtest (X, Y, NREPS, FUNC)
% -- Function File: PVAL = randtest (X, Y, NREPS, FUNC, SEED)
% -- Function File: [PVAL, STAT] = randtest (...)
% -- Function File: [PVAL, STAT, FPR] = randtest (...)
% -- Function File: [PVAL, STAT, FPR, PERMSTAT] = randtest (...)
%
%     'PVAL = randtest (X, Y)' uses the approach of Manly [1] to perform
%     a randomization (or permutation) test of the null hypothesis that
%     coefficients from the regression of Y on X are significantly different
%     from 0. The value returned is a 2-tailed p-value. Note that the Y values
%     are centered before randomization or permutation to also provide valid
%     null hypothesis tests of the intercept. To include an intercept term in
%     the regression, X must contain a column of ones.
%
%     Hint: For one-sample or two-sample randomization/permutation tests,
%     please use the 'randtest1' or 'randtest2' functions respectively.
%
%     'PVAL = randtest (X, Y, NREPS)' specifies the number of resamples without
%     replacement to take in the randomization test. By default, NREPS is 5000.
%     If the number of possible permutations is smaller than NREPS, the test
%     becomes exact. For example, if the number of sampling units (i.e. rows
%     in Y) is 6, then the number of possible permutations is factorial (6) =
%     720, so NREPS will be truncated at 720 and sampling will systematically
%     evaluate all possible permutations. 
%
%     'PVAL = randtest (X, Y, NREPS, FUNC)' also specifies a custom function
%     calculated on the original samples, and the permuted or randomized
%     resamples. Note that FUNC must compute statistics related to regression,
%     and should either be a:
%        o function handle or anonymous function,
%        o string of function name, or
%        o a cell array where the first cell is one of the above function
%          definitions and the remaining cells are (additional) input arguments 
%          to that function (other than the data arguments).
%        See the built-in demos for example usage with @mldivide for linear
%        regression coefficients, or with @cor for the correlation coefficient.
%        The default value of FUNC is @mldivide.
%
%     'PVAL = randtest (X, Y, NREPS, FUNC, SEED)' initialises the Mersenne
%     Twister random number generator using an integer SEED value so that
%     the results of 'randtest' results are reproducible when the
%     test is approximate (i.e. when using randomization if not all permutations
%     can be evaluated systematically).
%
%     '[PVAL, STAT] = randtest (...)' also returns the test statistic.
%
%     '[PVAL, STAT, FPR] = randtest (...)' also returns the minimum false
%     positive risk (FPR) calculated for the p-value, computed using the
%     Sellke-Berger approach.
%
%     '[PVAL, STAT, FPR, PERMSTAT] = randtest (...)' also returns the
%     statistics of the permutation distribution.
%
%  Bibliography:
%  [1] Manly (1997) Randomization, Bootstrap and Monte Carlo Method in Biology.
%       2nd Edition. London: Chapman & Hall.
%
%  randtest (version 2024.04.17)
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
%  along with this program.  If not, see http://www.gnu.org/licenses/

function [pval, stat, fpr, STATS] = randtest (x, y, nreps, func, seed)


  % Check if we are running Octave or Matlab
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

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
  if (nargin < 2)
    error ('randtest: X and Y must be provided')
  end
  if (nargin > 6)
    error ('randtest: Too many input arguments')
  end
  if (nargout > 4)
    error ('randtest: Too many output arguments')
  end
  
  % Set defaults
  if ( (nargin < 3) || isempty (nreps) )
    nreps = 5000;
  end
  if ( (nargin > 3) && (~ isempty (func)) )
    if (iscell (func))
      args = func(2:end);
      if (ischar (func{1}))
        % Convert character string of a function name to a function handle
        func = str2func (func{1});
      else
        func = func{1};
      end
      func = @(varargin) func (varargin{:}, args{:});
    elseif (ischar (func))
      % Convert character string of a function name to a function handle
      func = str2func (func);
    elseif (isa (func, 'function_handle'))
      % Do nothing
    else
      error ('randtest: FUNC must be a function name or function handle')
    end
  else
    func = @mldivide;
  end
  if ( (nargin > 4) && (~ isempty (seed)) )
    % Set random seed
    rand ('seed', seed);
  end

  % Remove NaN data values
  ridx = any (isnan (cat (2, x, y)), 2);
  x(ridx, :) = [];
  y(ridx, :) = [];

  % Get size of the data
  szx = size (x);
  szy = size (y);
  if (numel (szx) > 2)
    error ('randtest: X has too many dimensions')
  end
  if (numel (szy) > 2)
    error ('randtest: Y has too many dimensions')
  end
  if (szy(2) > 1)
    error ('randtest: Y must be a column vector')
  end
  if (szx(1) ~= szy(1))
    error (cat (2, 'randtest: X must be a vector or matrix with the same', ... 
                    'number of rows as y'))
  end
  n = szy (1);

  % Check for infinite values in the data
  if ( any (isinf (x(:))) || any (isinf (y)) )
    error ('randtest: X and Y cannot not contain inf values')
  end

  % Compute test statistic on the original data
  stat = func (x, y);

  % Centre y-values (so that intercept under the null is 0)
  y = y - mean (y);

  % Evaluate return value of FUNC
  sz_stat = size (stat);
  p = sz_stat(1);
  if (sz_stat(2) > 1)
    error ('randtest: FUNC (X, Y) must return a scalar or column vector')
  end

  % Create permutations or randomized samples
  if (factorial (n) < nreps)
    I = perms (1:n).';                  % For exact (permutation) test
    nreps = factorial (n);
  else 
    [jnk, I] = sort (rand (n, nreps));  % For approximate (randomization) test
  end
  Y = y(I);

  % Check if we can vectorize function evaluations
  VECTORIZED = all (bsxfun (@eq, size (func (x, repmat (y, 1, 2))), [p, 2]));

  % Perform function evaluations
  if VECTORIZED
    STATS = func (x, Y);
  else
    if (PARALLEL)
      if (ISOCTAVE)
        STATS = parcellfun (inf, func, num2cell (repmat(x, 1, nreps), 1), ...
                             num2cell (Y, 1));
      else
        STATS = zeros (1, nreps);
        parfor b = 1:nreps
          STATS(b) = func (x, Y(:, b))
        end
      end
    else
      STATS = cellfun (func, num2cell (repmat(x, 1, nreps), 1), ...
                             num2cell (Y, 1));
    end
  end

  % Calculate two-tailed p-value(s) by linear interpolation
  pval = nan (p, 1);
  res_lim = 1 / nreps;
  for j = 1:p
    if (~ isnan (stat(j)))
      [u, jnk, P] = bootcdf (abs (STATS(j,:)), true);
      if (numel (u) > 1)
        if (abs (stat(j)) < u(1))
          pval(j) = interp1 (u, P, abs (stat(j)), 'linear', 1);
        else
          pval(j) = interp1 (u, P, abs (stat(j)), 'linear', res_lim);
        end
      else 
        pval(j) = P;
      end
    end
  end

  % Compute minimum false positive risk
  if (nargout > 2)
    fpr = pval2fpr (pval);
  end

end

%--------------------------------------------------------------------------

% FUNCTION TO COMPUTE MINIMUM FALSE POSITIVE RISK (FPR)

function fpr = pval2fpr (p)

  % Subfunction to compute minimum false positive risk. These are calculated
  % from a Bayes factor based on the sampling distributions of the p-value and
  % that H0 and H1 have equal prior probabilities. This is called the Sellke-
  % Berger approach.
  % 
  % References:
  %  Held and Ott (2018) On p-Values and Bayes Factors. 
  %    Annu. Rev. of Stat. Appl. 5:393-419
  %  David Colquhoun (2019) The False Positive Risk: A Proposal 
  %    Concerning What to Do About p-Values, The American Statistician, 
  %    73:sup1, 192-201, DOI: 10.1080/00031305.2018.1529622 

  % Calculate minimum Bayes Factor (P(H0) / P(H1)) by the Sellke-Berger method 
  logp = min (log (p), -1);
  minBF = exp (1 + logp + log (-logp));

  % Calculate the false-positive risk from the minumum Bayes Factor
  L10 = 1 ./ minBF;      % Convert to Maximum Likelihood ratio L10 (P(H1)/P(H0))
  fpr = max (0, 1 ./ (1 + L10));  % Calculate minimum false positive risk 
  fpr(isnan(p)) = NaN; 

end

%--------------------------------------------------------------------------

%!demo
%!
%! % Randomization or permutation test for linear regression (without intercept)
%! % cd4 data in DiCiccio and Efron (1996) Statistical Science
%! X = [212 435 339 251 404 510 377 335 410 335 ...
%!      415 356 339 188 256 296 249 303 266 300]';
%! Y = [247 461 526 302 636 593 393 409 488 381 ...
%!      474 329 555 282 423 323 256 431 437 240]';
%!
%! % Randomization test to assess the statistical significance of the
%! % regression coefficient being different from 0. (Model: y ~ x or y = 0 + x,
%! % i.e. linear regression through the origin)
%! [pval, stat]  = randtest (X, Y, 5000) % Default value of FUNC is @mldivide
%!

%!demo
%!
%! % Randomization or permutation test for linear regression (with intercept)
%! % cd4 data in DiCiccio and Efron (1996) Statistical Science
%! X = [212 435 339 251 404 510 377 335 410 335 ...
%!      415 356 339 188 256 296 249 303 266 300]';
%! Y = [247 461 526 302 636 593 393 409 488 381 ...
%!      474 329 555 282 423 323 256 431 437 240]';
%! N = numel (Y);
%!
%! % Randomization test to assess the statistical significance of the
%! % regression coefficients (intercept and slope) being different from 0.
%! % (Model: y ~ 1 + x, i.e. linear regression with intercept)
%! X1 = cat (2, ones (N, 1), X);
%! [pval, stat]  = randtest (X1, Y, 5000) % Default value of FUNC is @mldivide

%!demo
%!
%! % Randomization or permutation test for the correlation coefficient
%! % cd4 data in DiCiccio and Efron (1996) Statistical Science
%! X = [212 435 339 251 404 510 377 335 410 335 ...
%!      415 356 339 188 256 296 249 303 266 300]';
%! Y = [247 461 526 302 636 593 393 409 488 381 ...
%!      474 329 555 282 423 323 256 431 437 240]';
%!
%! % Randomization test to assess the statistical significance of the
%! % correlation coefficient being different from 0. This is equivalent to 
%! % the slope regression coefficient for a linear regression (with intercept)
%! % of standardized x and y values.
%! [pval, stat] = randtest (X, Y, 5000, @cor)

%!test
%!
%! % Test various capabilities of randtest
%! X = randn (3,1);
%! Y = randn (3,1);
%! pval1 = randtest (X, Y);
%! pval2 = randtest (X, Y, 500);
%! randtest (X, Y, [], []);
%! X = randn (9,1);
%! Y = randn (9,1);
%! pval3 = randtest (X, Y, 5000);
%! pval4 = randtest (X, Y, [], [], 1);
%! pval5 = randtest (X, Y, [], @mldivide, 1);
%! pval6 = randtest (X, Y, [], @cor, 1);
