% Computes percentile confidence interval(s) directly from a vector (or row-
% major matrix) of bootstrap statistics.
%
% -- Function File: CI = bootint (BOOTSTAT)
% -- Function File: CI = bootint (BOOTSTAT, PROB)
% -- Function File: CI = bootint (BOOTSTAT, PROB, ORIGINAL)
%
%     'CI = bootint (BOOTSTAT)' computes simple 95% percentile confidence
%     intervals [1,2] directly from the vector, or rows* of the matrix in
%     BOOTSTAT, where BOOTSTAT contains bootstrap statistics such as those
%     generated using the `bootstrp` function. Depending on the application,
%     bootstrap confidence intervals with better coverage and accuracy can
%     be computed using the various dedicated bootstrap confidence interval
%     functions from the statistics-resampling package.
%
%        * The matrix should have dimensions P * NBOOT, where P corresponds to
%          the number of parameter estimates and NBOOT corresponds to the number
%          of bootstrap samples.
%
%     'CI = bootint (BOOTSTAT, PROB)' returns confidence intervals, where
%     PROB is numeric and sets the lower and upper bounds of the confidence
%     interval(s). The value(s) of PROB must be between 0 and 1. PROB can
%     either be:
%       <> scalar: To set the central mass of normal confidence intervals
%                  to 100*PROB%
%       <> vector: A pair of probabilities defining the lower and upper
%                  percentiles of the confidence interval(s) as 100*(PROB(1))%
%                  and 100*(PROB(2))% respectively.
%          The default value of PROB is the vector: [0.025, 0.975], for an
%          equal-tailed 95% percentile confidence interval.
%
%     'CI = bootint (BOOTSTAT, PROB, ORIGINAL)' uses the ORIGINAL estimates
%     associated with BOOTSTAT to correct PROB and the resulting confidence
%     intervals (CI) for median bias. The confidence intervals returned in CI
%     therefore become bias-corrected percentile intervals [3,4].
%
%  BIBLIOGRAPHY:
%  [1] Efron (1979) Bootstrap Methods: Another look at the jackknife.
%        Annals Stat. 7,1-26
%  [2] Efron, and Tibshirani (1993) An Introduction to the Bootstrap. 
%        New York, NY: Chapman & Hall
%  [3] Efron (1981) Nonparametric Standard Errors and Confidence Intervals.
%        Can J Stat. 9(2):139-172
%  [4] Efron (1982) The jackknife, the bootstrap, and other resampling plans.
%        SIAM-NSF, CBMS #38
%
%  bootint (version 2024.05.19)
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

function CI = bootint (Y, PROB, T0)

  % Check input and output arguments
  if (nargin < 2)
    PROB = 0.95;
  end
  if (nargin > 3)
    error ('bootint: Too many input arguments.')
  end
  if (nargout > 1)
    error ('bootint: Too many output arguments.')
  end

  % Evaluate the dimensions of y
  sz = size (Y);
  if all (sz == 1)
    error ('bootint: Y must be either a vector or a P * NBOOT matrix')
  end
  if (sz(2) == 1)
    p = 1;
    Y = Y.';
    nboot = numel (Y);
  else
    p = sz(1);
    nboot = sz(2);
    if (p > nboot)
      warning ('bootint: The dimensions of the matrix in Y should be P * NBOOT')
    end
  end

  % Evaluate PROB input argument
  if ( (nargin < 2) || (isempty (PROB)) )
    PROB = [0.025, 0.975];
  else
    nprob = numel (PROB);
    if (~ isa (PROB, 'numeric') || (nprob > 2))
      error ('bootint: PROB must be a scalar or a vector of length 2')
    end
    if (size (PROB, 1) > 1)
      PROB = PROB.';
    end
    if (any ((PROB < 0) | (PROB > 1)))
      error ('bootint: Value(s) in PROB must be between 0 and 1')
    end
    if (nprob > 1)
      % PROB is a pair of probabilities
      % Make sure probabilities are in the correct order
      if (PROB(1) > PROB(2) )
        error (cat (2, 'bootint: The pair of probabilities must be in', ...
                       ' ascending numeric order'))
      end
    else
      PROB = arrayfun (@(c) (1 + c * PROB) / 2, [-1, 1]);
    end
  end

  % Apply median bias correction to PROB if T0 is provided
  if (nargin > 2)
    % Error checking for T0
    if (~ isa (PROB, 'numeric'))
      error ('bootint: ORIGINAL must be numeric')
    end
    szT0 = size (T0);
    if ( (all (szT0 > 1)) || (szT0(2) > 1) )
      error ('bootint: ORIGINAL must be a scalar or column vector')
    end
    if (p > 1)
      if (szT0(1) ~= p)
        error ('bootint: BOOTSTAT and ORIGINAL must have the same number of rows')
      end
    else
      if (szT0(1) > 1)
        error ('bootint: If BOOTSTAT is a vector, ORIGINAL must be a scalar')
      end
    end
    % Create distribution functions
    stdnormcdf = @(x) 0.5 * (1 + erf (x / sqrt (2)));
    stdnorminv = @(p) sqrt (2) * erfinv (2 * p - 1);
    % Calculate the median bias correction constant (z0)
    z0 = stdnorminv (sum (bsxfun (@lt, Y, T0), 2) / nboot);
    if (~ all (isfinite (z0)))
      % Revert to percentile bootstrap confidence intervals
      warning ('bootint:biasfail', ...
               cat (2, 'Unable to calculate the bias correction constant;', ...
                       ' reverting to simple percentile intervals.\n'))
      z0 = zeros (p, 1);
    end
    % Calculate BC percentiles
    z = stdnorminv (PROB);
    PROB = stdnormcdf (bsxfun (@plus, 2 * z0, z));
  end

  % Compute confidence intervals
  CI = nan (p, 2);
  for j = 1:p
    [t1, cdf] = bootcdf (Y(j, :), true, 1);
    CI(j, 1) = interp1 (cdf, t1, PROB(1) , 'linear', min (t1));
    CI(j, 2) = interp1 (cdf, t1, PROB(2) , 'linear', max (t1));
  end
  CI(:, isnan(PROB)) = NaN;
  if (PROB(1) == 0)
    CI(:, 1) = -inf;
  end
  if (PROB(2) == 1)
    CI(:, 2) = +inf;
  end

end

%!demo
%!
%! % Law school data
%! data = [576, 3.39; 635, 3.30; 558, 2.81; 578, 3.03; 666, 3.44; ...
%!         580, 3.07; 555, 3.00; 661, 3.43; 661, 3.36; 605, 3.13; ...
%!         653, 3.12; 575, 2.74; 545, 2.76; 572, 2.88; 594, 2.96];
%! x = data(:, 1);
%! y = data(:, 2);
%! r = cor (x, y);
%!
%! % 95% confidence interval for the mean 
%! bootstat = bootstrp (4999, @cor, x, y);
%! CI_per  = bootint (bootstat,0.95)    % 95% simple percentile interval
%! CI_cper = bootint (bootstat,0.95,r)  % 95% bias-corrected percentile interval
%!
%! % Please be patient, the calculations will be completed soon...

%!test
%!
%! % Law school data
%! data = [576, 3.39; 635, 3.30; 558, 2.81; 578, 3.03; 666, 3.44; ...
%!         580, 3.07; 555, 3.00; 661, 3.43; 661, 3.36; 605, 3.13; ...
%!         653, 3.12; 575, 2.74; 545, 2.76; 572, 2.88; 594, 2.96];
%! x = data(:, 1);
%! y = data(:, 2);
%! r = cor (x, y);
%!
%! % 95% confidence interval for the mean 
%! bootstat = bootstrp (4999, @cor, x, y);
%! CI = bootint (bootstat,0.95);             % 95% percentile interval
%! CI = bootint (bootstat,[0.025,0.975]);    % 95% percentile interval
%! CI = bootint (bootstat,0.95,r);           % 95% bias-corrected interval
%! CI = bootint (bootstat,[0.025,0.975],r);  % 95% bias-corrected interval
