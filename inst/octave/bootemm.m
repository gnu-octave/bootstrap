%  Function File: bootemm
%
%  EMMEANS = bootemm (STATS, DIM)
%  EMMEANS = bootemm (STATS, DIM, NBOOT)
%  EMMEANS = bootemm (STATS, DIM, NBOOT, ALPHA)
%  EMMEANS = bootemm (STATS, DIM, NBOOT, ALPHA, NPROC)
%  EMMEANS = bootemm (STATS, DIM, NBOOT, ALPHA, NPROC, SEED)
%  bootemm (STATS, DIM)
%  bootemm (STATS, DIM, ...)
%
%  Non-parametric bootstrap of the estimated marginal means from a linear model.
%  bootemm accepts as input the STATS structure from fitlm or anovan functions
%  (from the v1.5+ of the Statistics package in OCTAVE) and returns a structure,
%  EMMEANS, which contains the following fields:
%    original: contains the estimated marginal means from the original model
%    bias: contains the bootstrap estimate of bias
%    std_error: contains the bootstrap standard error
%    CI_lower: contains the lower bound of the 95% bootstrap confidence interval
%    CI_upper: contains the upper bound of the 95% bootstrap confidence interval
%  The method uses bootknife resampling [1], which involves creating leave-one-
%  out jackknife samples of size n - 1 from the n residuals and then drawing
%  samples of size n with replacement from the jackknife samples. The resampling
%  of residuals is also balanced in order to reduce bias and Monte Carlo error
%  [2,3]. By default, the confidence intervals constructed are bias-corrected
%  percentile intervals [4]. The list of means and their bootstrap statistics
%  correspond to the names STATS.grpnames. If no output is requested, the
%  results are printed to stdout.
%
%  EMMEANS = bootemm (STATS, NBOOT) also specifies the number of bootstrap 
%  samples. NBOOT can be a scalar value (for single bootstrap), or vector of
%  upto two positive integers (for double bootstrap). By default, NBOOT is
%  [2000,0]. See documentation for the bootknife function for more information.
%
%  EMMEANS = bootemm (STATS, NBOOT, ALPHA) where ALPHA is numeric and
%  sets the lower and upper bounds of the confidence interval(s). The value(s)
%  of ALPHA must be between 0 and 1. ALPHA can either be a scalar value to set
%  the (nominal) central coverage, or a vector of 2 numeric values corresponding
%  to a pair of probabilities to set the (nominal) lower and upper percentiles.
%  The value(s) of ALPHA must be between 0 and 1. The method for constructing
%  confidence intervals is determined by the combined settings of ALPHA and
%  NBOOT. See documentation for the bootknife function for more information. 
%  Confidence intervals are not calculated when the value(s) of ALPHA is/are NaN. 
%
%  EMMEANS = bootemm (STATS, NBOOT, ALPHA, NPROC) also sets the number of
%  parallel processes to use to accelerate computations on multicore machines.
%  This feature requires the Parallel package (in Octave). See documentation
%  for the bootknife function for more information.
%
%  EMMEANS = bootemm (STATS, NBOOT, ALPHA, NPROC, SEED) also sets the random
%  SEED for the random number generator used for the resampling. This feature
%  can be used to make the results of the bootstrap reproducible.
%
%  Requirements: The function file boot.m (or better boot.mex) and bootknife
%  also distributed in the statistics-bootstrap package. bootcoeff is only
%  supported in GNU Octave and requires the Statistics package version 1.5+.
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
%
%  bootemm (version 2023.01.12)
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


function emmeans = bootemm (STATS, dim, nboot, alpha, ncpus, seed)

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
    alpha = 0.05;
  end
  if (nargin < 5)
    ncpus = 0;
  end
  if (nargin > 4)
    boot (1, 1, false, seed);
  end

  % Error checking
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));
  if (~ ISOCTAVE)
    error ('bootemm: only supported by Octave')
  end
  statspackage = ismember ({info.Name}, 'statistics');
  if (~ any (statspackage)) || (str2num (info (statspackage).Version(1:3)) < 1.5)
    error ('bootemm: requires version >= 1.5 of the statistics package')
  end
  if (ismember (dim, find (STATS.continuous)))
    error ('bootemm: estimated marginal means are only calculated for categorical variables')
  end

  % Fetch required information from STATS structure
  X = full (STATS.X);
  b = STATS.coeffs(:,1);
  fitted = X * b;
  W = full (STATS.W);
  w = diag (W);
  se = diag (W).^(-0.5);
  resid = STATS.resid;   % weighted residuals
  y = fitted + resid .* se;
  n = numel (resid);

  % Prepare the hypothesis matrix (H)
  Nd = numel (dim);
  df = STATS.df;
  i = 1 + cumsum(df);
  k = find (sum (STATS.terms(:,dim), 2) == sum (STATS.terms, 2));
  Nb = 1 + sum(df(k));
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

  % Define bootfun and data for case resampling of raw data
  % Robust to violations of homoskedasticity and normality assumptions
  bootfun = @(X, y, w) H * pinv (diag (w) * X) * (w .* y);
  data = {X, y, w};

  % Perform bootstrap
  warning ('off','bootknife:lastwarn')
  if nargout > 0
    emmeans = bootknife (data, nboot, bootfun, alpha, [], ncpus);
  else
    bootknife (data, nboot, bootfun, alpha, [], ncpus);
  end
  warning ('on','bootknife:lastwarn')

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