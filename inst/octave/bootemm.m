% -- Function File: bootemm (STATS, DIM)
% -- Function File: bootemm (STATS, DIM, NBOOT)
% -- Function File: bootemm (STATS, DIM, NBOOT, ALPHA)
% -- Function File: bootemm (STATS, DIM, NBOOT, ALPHA)
% -- Function File: bootemm (STATS, DIM, NBOOT, ALPHA, SEED)
% -- Function File: CI = bootemm (STATS, DIM, ...)
% -- Function File: [CI, BOOTSTAT] = bootemm (STATS, DIM, ...)
%
%     'bootemm (STATS, DIM)' uses the STATS structure output from fitlm or
%     anovan functions (from the v1.5+ of the Statistics package in OCTAVE)
%     and Bayesian bootstrap to compute and return the following statistics
%     along the dimension DIM:
%        • original: estimated marginal means (EMMs) from the original data
%        • bias: bootstrap estimate of the bias of the EMMs
%        • std_error: bootstrap estimate of the standard error of the EMMs
%        • CI_lower: lower bound of the 95% bootstrap confidence interval
%        • CI_upper: upper bound of the 95% bootstrap confidence interval
%          Here, the confidence intervals, or credible intervals in the context
%          of the Bayesian statistical framework, are percentile intervals [2].
%
%     'bootemm (STATS, DIM, NBOOT)' specifies the number of bootstrap resamples,
%     where NBOOT must be a positive integer. If empty, tHe default value of
%     NBOOT is the scalar: 2000.
%
%     'bootemm (STATS, DIM, NBOOT, ALPHA)', where ALPHA is numeric and sets the
%     the lower and upper bounds of the confidence interval(s). The value(s) of
%     ALPHA must be between 0 and 1. ALPHA can either be:
%        • scalar: To set the (nominal) central coverage of equal-tailed
%                  percentile confidence intervals to 100*(1-ALPHA)%.
%        • vector: A pair of probabilities defining the (nominal) lower and
%                  upper percentiles of the confidence interval(s) as
%                  100*(ALPHA(1))% and 100*(ALPHA(2))% respectively. 
%        Confidence intervals are not calculated when the value(s) of ALPHA
%        is/are NaN. The default value of  ALPHA is the vector: [.025, .975], 
%        for a 95% confidence interval.
%
%     'bootemm (STATS, DIM, NBOOT, ALPHA, SEED)', initialises the Mersenne
%     Twister random number generator using an integer SEED value so that
%     bootemm results are reproducible.
%
%     'CI = bootemm (STATS, DIM, ...) returns the confidence intervals of the
%     estimated marginal means for the linear model.
%
%     '[CI, BOOTSTAT] = bootemm (STATS, DIM, ...) also returns the bootstrap
%     statistics for the estimated marginal means.
%
%  Requirements: bootemm is only supported in GNU Octave and requires the
%  Statistics package version 1.5+.
%
%  Bibliography:
%  [1] Rubin (1981) The Bayesian Bootstrap. Ann. Statist. 9(1):130-134
%  [2] Efron and Tibshirani. Chapter 16 Hypothesis testing with the
%       bootstrap in An introduction to the bootstrap (CRC Press, 1994)
%
%  bootemm (version 2023.05.02)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/
%
%  Copyright 2019 Andrew Charles Penn
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of  the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.


function [emm, bootstat] = bootemm (STATS, dim, nboot, alpha, seed)

  % Check input aruments
  if (nargin < 2)
    error ('bootemm usage: ''bootemm (STATS, dim)'' atleast 2 input arguments required');
  end
  if (nargin < 3)
    nboot = 2000;
  end
  if (nargin < 3)
    nboot = 2000;
  end
  if (nargin < 4)
    alpha = [0.025, 0.975];
  end
  if (nargin > 4)
    if (ISOCTAVE)
      rande ('seed', seed);
    else
      rng ('default');
    end
  end

  % Error checking
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));
  if (~ ISOCTAVE)
    error ('bootemm: Only supported by Octave')
  end
  statspackage = ismember ({info.Name}, 'statistics');
  if ((~ any (statspackage)) || (str2double (info (statspackage).Version(1:3)) < 1.5))
    error ('bootemm: Requires version >= 1.5 of the statistics package')
  end
  if (ismember (dim, find (STATS.continuous)))
    error ('bootemm: Estimated marginal means are only calculated for categorical variables')
  end

  % Fetch required information from STATS structure
  X = full (STATS.X);
  b = STATS.coeffs(:,1);
  fitted = X * b;
  resid = STATS.resid;
  y = fitted + resid;
  n = numel (resid);

  % Prepare the hypothesis matrix (H)
  df = STATS.df;
  i = 1 + cumsum(df);
  k = find (sum (STATS.terms(:,dim), 2) == sum (STATS.terms, 2));
  Nt = numel (k);
  L = zeros (n, sum (df) + 1);
  for j = 1:Nt
    L(:, i(k(j)) - df(k(j)) + 1 : i(k(j))) = STATS.X(:,i(k(j)) - ...
                                             df(k(j)) + 1 : i(k(j)));
  end
  L(:,1) = 1;
  H = unique (L, 'rows', 'stable');
  Ng = size (H, 1);
  idx = zeros (Ng, 1);
  for k = 1:Ng
    idx(k) = find (all (L == H(k, :), 2),1);
  end

  % Perform Bayesian bootstrap
  switch (nargout)
    case 0
      bootbayes (y, X, nboot, alpha, [], H);
    case 1
      emm = bootbayes (y, X, nboot, alpha, [], H);
    otherwise
      [emm, bootstat] = bootbayes (y, X, nboot, alpha, [], H);
  end

end

%!demo
%!
%! dv =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! g = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%!
%! [P,ATAB,STATS] = anovan (dv,g);
%! DIM = 1;
%! STATS.grpnames{DIM}
%! bootemm (STATS, DIM);              # 95% confidence intervals 

%!demo
%!
%! dv =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! g = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%!
%! [P,ATAB,STATS] = anovan (dv,g);
%! DIM = 1;
%! STATS.grpnames{DIM}
%! bootemm (STATS, DIM, 2000, 0.166); # 83.4% confidence intervals (overlap at p > 0.05)