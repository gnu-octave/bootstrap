% Computes percentile confidence interval(s) from a vector (or row-major
% matrix) of bootstrap statistics.
%
% -- Function File: CI = bootint (Y)
% -- Function File: CI = bootint (Y, PROB)
%
%     'CI = bootint (Y)' computes 95% percentile confidence intervals from
%     the vector, or rows* of the matrix, Y, where Y contains bootstrap
%     statistics, such as those generated using the `bootstrp` function.
%     Depending on the application, bootstrap confidence intervals with better
%     coverage and accuracy can be computed using the various dedicated
%     bootstrap functions from the statistics-resampling package.
%
%        * The matrix should have dimensions P * NBOOT, where P corresponds to
%          the number of parameter estimates and NBOOT corresponds to the number
%          of bootstrap samples.
%
%     'CI = bootint (Y, PROB)' returns confidence intervals, where PROB is
%     numeric and sets the lower and upper bounds of the confidence interval(s).
%     The value(s) of PROB must be between 0 and 1. PROB can either be:
%       <> scalar: To set the central mass of normal confidence intervals
%                  to 100*PROB%
%       <> vector: A pair of probabilities defining the lower and upper
%                  percentiles of the confidence interval(s) as 100*(PROB(1))%
%                  and 100*(PROB(2))% respectively.
%          The default value of PROB is the vector: [0.025, 0.975], for an
%          equal-tailed 95% percentile confidence interval.
%
%  bootint (version 2024.04.24)
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

function CI = bootint (Y, PROB)

  % Check input and output arguments
  if (nargin < 2)
    PROB = 0.95;
  end
  if (nargin > 2)
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

%!test
%!
%! % Simulate (log-normal) data
%! randn ('seed', 1);
%! Y = exp (randn (5, 999));
%!
%! % 95% confidence interval for the mean 
%! CI = bootint (Y,0.95);          # 95% percentile interval
%! CI = bootint (Y,[0.025,0.975]); # 95% percentile interval
