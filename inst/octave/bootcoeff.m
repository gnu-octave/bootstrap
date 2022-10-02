%  Function File: bootcoeff
%
%  coeffs = bootcoeff (stats)
%  coeffs = bootcoeff (stats, nboot)
%  coeffs = bootcoeff (stats, nboot, alpha)
%  coeffs = bootcoeff (stats, nboot, alpha, nproc)
%  coeffs = bootcoeff (stats, nboot, alpha, nproc, seed)
%
%  Semi-parametric bootstrap of the regression coefficients from a linear model.
%  bootcoeff accepts as input the stats structure from fitlm or anovan functions
%  (from the v1.5+ of the Statistics package in OCTAVE) and returns a structure,
%  coeffs, which contains the following fields:
%    original: contains the regression coefficients from the original data
%    bias: contains the bootstrap estimate of bias
%    std_error: contains the bootstrap standard error
%    CI_lower: contains the lower bound of the bootstrap confidence interval
%    CI_upper: contains the upper bound of the bootstrap confidence interval
%  By default, the confidence intervals are 95% bias-corrected intervals. If
%  no output is requested, the results are printed to stdout. The list of
%  coefficients and their bootstrap statistics correspond to the names in
%  stats.coeffnames, which are defined by the contrast coding in stats.contrasts.
%  The rows of stats.contrasts correspond to the names in stats.grpnames.
%
%  coeffs = bootcoeff (stats, nboot) also specifies the number of bootstrap
%  samples. nboot must be a scalar. By default, nboot is 2000.
%
%  coeffs = bootcoeff (stats, nboot, alpha) also sets the lower and upper
%  confidence interval ends. The value(s) in alpha must be between 0 and 1.
%  If alpha is a scalar value, the nominal lower and upper percentiles of
%  the confidence are 100*(alpha/2)% and 100*(1-alpha/2)% respectively, and
%  the intervals are bias-corrected with nominal central coverage 100*(1-alpha)%.
%  If alpha is a vector with two elements, alpha becomes the quantiles for
%  percentile bootstrap confidence intervals. If alpha is empty, NaN is returned
%  for the confidence interval ends. The default value for alpha is 0.05. 
%
%  coeffs = bootcoeff (stats, nboot, alpha, nproc) also sets the number of
%  parallel processes to use to accelerate computations on multicore machines.
%  This feature requires the Parallel package (in Octave).
%
%  coeffs = bootcoeff (stats, nboot, alpha, nproc, seed) also sets the random
%  seed for the random number generator used for the resampling. This feature
%  can be used to make the results of the bootstrap reproducible.
%
%  bootcoef is only supported in GNU Octave and requires the Statistics package
%  version 1.5 or later.
%
%  bootcoeff v0.1.1.0 (30/09/2022)
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
    boot (1, 1, true, [], seed);
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
    [coeffs, bootstat] = bootknife (resid, nboot, bootfun, alpha, [], ncpus);
    warning ('on','bootknife:lastwarn')
  else
    coeffs = [];
    bootknife (resid, nboot, bootfun, alpha, [], ncpus);
  end

end