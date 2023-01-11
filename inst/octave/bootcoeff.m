%  Function File: bootcoeff
%
%  COEFFS = bootcoeff (STATS)
%  COEFFS = bootcoeff (STATS, NBOOT)
%  COEFFS = bootcoeff (STATS, NBOOT, ALPHA)
%  COEFFS = bootcoeff (STATS, NBOOT, ALPHA, NPROC)
%  COEFFS = bootcoeff (STATS, NBOOT, ALPHA, NPROC, SEED)
%  bootcoeff (STATS)
%  bootcoeff (STATS, ...)
%
%  Semi-parametric bootstrap of the regression coefficients from a linear model.
%  bootcoeff accepts as input the STATS structure from fitlm or anovan functions
%  (from the v1.5+ of the Statistics package in OCTAVE) and returns a structure,
%  COEFFS, which contains the following fields:
%    original: contains the regression coefficients from the original data
%    bias: contains the bootstrap estimate of bias
%    std_error: contains the bootstrap standard error
%    CI_lower: contains the lower bound of the 95% bootstrap confidence interval
%    CI_upper: contains the upper bound of the 95% bootstrap confidence interval
%  The method uses bootknife resampling [1], which involves creating leave-one-
%  out jackknife samples of size n - 1 from the n residuals and then drawing
%  samples of size n with replacement from the jackknife samples. The resampling
%  of residuals is also balanced in order to reduce bias and Monte Carlo error
%  [2,3]. By default, the confidence intervals constructed are bias-corrected
%  percentile intervals [4,5]. The list of coefficients and their bootstrap
%  statistics correspond to the names in STATS.coeffnames, which are defined by
%  the contrast coding in STATS.contrasts. The rows of STATS.contrasts 
%  correspond to the names in STATS.grpnames. If no output is requested, the
%  results are printed to stdout.
%  
%  COEFFS = bootcoeff (STATS, NBOOT) also specifies the number of bootstrap
%  samples. NBOOT must be a scalar. By default, NBOOT is 2000.
%
%  COEFFS = bootcoeff (STATS, NBOOT, ALPHA) where ALPHA is numeric and
%  sets the lower and upper bounds of the confidence interval(s). The value(s)
%  of ALPHA must be between 0 and 1. ALPHA can either be:
%
%  1) a scalar value to set the (nominal) central coverage to 100*(1-ALPHA)%
%  with (nominal) lower and upper percentiles of the confidence intervals at
%  100*(ALPHA/2)% and 100*(1-ALPHA/2)% respectively. The intervals constructed
%  bias-corrected [4].
%
%  2) a vector containing a pair of probabilities to set the (nominal) lower and
%  upper percentiles of the confidence interval(s) at 100*(ALPHA(1))% and
%  100*(ALPHA(2))%. The intervals constructed are simple percentile intervals.
%
%  Confidence interval endpoints are not calculated when the value(s) of ALPHA
%  is/are NaN. If empty (or not specified), the default value for ALPHA is 0.05
%
%  COEFFS = bootcoeff (STATS, NBOOT, ALPHA, NPROC) also sets the number of
%  parallel processes to use to accelerate computations on multicore machines.
%  This feature requires the Parallel package (in Octave).
%
%  COEFFS = bootcoeff (STATS, NBOOT, ALPHA, NPROC, SEED) also sets the random
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
%  [4] Efron (1981) Nonparametric Standard Errors and Confidence Intervals.
%        Can J Stat. 9(2:139-158
%
%  bootcoeff (version 2022.12.16)
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


function coeffs = bootcoeff (stats, nboot, alpha, ncpus, seed)

  % Check input aruments
  if (nargin < 1)
    error ('bootcoeff usage: ''bootcoeff (stats)'' atleast 1 input arguments required');
  end
  if (nargin < 2)
    nboot = 2000;
  end
  if (nargin < 3)
    alpha = 0.05;
  end
  if (nargin < 4)
    ncpus = 0;
  end
  if (nargin > 4)
    boot (1, 1, false, seed);
  end

  % Error checking
  if numel(nboot) > 1
    error ('bootcoeff only supports single bootstrap resampling')
  end
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));
  if ~ISOCTAVE
    error ('bootcoeff is only supported by Octave')
  end
  statspackage = ismember ({info.Name}, 'statistics');
  if (~ any (statspackage)) || (str2num(info(statspackage).Version(1:3)) < 1.5)
    error ('bootcoeff requires version >= 1.5 of the statistics package')
  end

  % Fetch required information from stats structure
  X = stats.X;
  b = stats.coeffs(:,1);
  fitted = X * b;
  lmfit = stats.lmfit;
  W = full (stats.W);
  se = diag (W).^(-0.5);
  resid = stats.resid;   % weighted residuals

  % Define bootfun for bootstraping the model residuals and returning the regression coefficients
  bootfun = @(r) lmfit (X, fitted + r .* se, W);

  % Perform bootstrap
  if nargout > 0
    warning ('off','bootknife:lastwarn')
    coeffs = bootknife (resid, nboot, bootfun, alpha, [], ncpus);
    warning ('on','bootknife:lastwarn')
  else
    bootknife (resid, nboot, bootfun, alpha, [], ncpus);
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
%! [P,ATAB,STATS] = anovan (dv,g,'contrasts','simple');
%! STATS.coeffnames
%! bootcoeff (STATS)