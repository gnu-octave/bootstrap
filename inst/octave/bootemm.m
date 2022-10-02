%  Function File: bootemm
%
%  emmeans = bootemm (stats, dim)
%  emmeans = bootemm (stats, dim, nboot)
%  emmeans = bootemm (stats, dim, nboot, alpha)
%  emmeans = bootemm (stats, dim, nboot, alpha, nproc)
%  emmeans = bootemm (stats, dim, nboot, alpha, nproc, seed)
%
%  Semi-parametric bootstrap of the estimated marginal means from a linear model.
%  bootemm accepts as input the stats structure from fitlm or anovan functions
%  (from the v1.5+ of the Statistics package in OCTAVE) and returns a structure,
%  emmeans, which contains the following fields:
%    original: contains the estimated marginal means from the original model
%    bias: contains the bootstrap estimate of bias
%    std_error: contains the bootstrap standard error
%    CI_lower: contains the lower bound of the bootstrap confidence interval
%    CI_upper: contains the upper bound of the bootstrap confidence interval
%  By default, the confidence intervals are 95% bias-corrected intervals. If
%  no output is requested, the results are printed to stdout. The list of
%  means and their bootstrap statistics correspond to the names stats.grpnames.
%
%  emmeans = bootemm (stats, nboot) also specifies the number of bootstrap
%  samples. nboot must be a scalar. By default, nboot is 2000.
%
%  emmeans = bootemm (stats, nboot, alpha) also sets the lower and upper
%  confidence interval ends. The value(s) in alpha must be between 0 and 1.
%  If alpha is a scalar value, the nominal lower and upper percentiles of
%  the confidence are 100*(alpha/2)% and 100*(1-alpha/2)% respectively, and
%  the intervals are bias-corrected with nominal central coverage 100*(1-alpha)%.
%  If alpha is a vector with two elements, alpha becomes the quantiles for
%  percentile bootstrap confidence intervals. If alpha is empty, NaN is returned
%  for the confidence interval ends. The default value for alpha is 0.05. 
%
%  emmeans = bootemm (stats, nboot, alpha, nproc) also sets the number of
%  parallel processes to use to accelerate computations on multicore machines.
%  This feature requires the Parallel package (in Octave).
%
%  emmeans = bootemm (stats, nboot, alpha, nproc, seed) also sets the random
%  seed for the random number generator used for the resampling. This feature
%  can be used to make the results of the bootstrap reproducible.
%
%  bootemm is only supported in GNU Octave and requires the Statistics package
%  version 1.5 or later.
%
%  bootemm v0.1.1.0 (30/09/2022)
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


function emmeans = bootemm (stats, dim, nboot, alpha, ncpus, seed)

  % Check input aruments
  if (nargin < 2)
    error ('bootemm usage: ''bootemm (stats, dim)'' atleast 2 input arguments required');
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
    boot (1, 1, true, [], seed);
  end

  % Error checking
  if numel(nboot) > 1
    error ('bootemm only supports single bootstrap resampling')
  end
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));
  if ~ISOCTAVE
    error ('bootemm is only supported by Octave')
  end
  statspackage = ismember ({info.Name}, 'statistics');
  if (~ any (statspackage)) || (str2num(info(statspackage).Version(1:3)) < 1.5)
    error ('bootemm requires version >= 1.5 of the statistics package')
  end
  if (ismember (dim, find (stats.continuous)))
    error ('bootemm: estimated marginal means are only calculated for categorical variables')
  end
    
  % Fetch required information from stats structure
  X = stats.X;
  b = stats.coeffs(:,1);
  fitted = X * b;
  lmfit = stats.lmfit;
  W = full (stats.W);
  se = diag (W).^(-0.5);
  resid = stats.resid;   % weighted residuals

  % Prepare the hypothesis matrix (H)
  Nd = numel (dim);
  n = numel (stats.resid);
  df = stats.df;
  i = 1 + cumsum(df);
  k = find (sum (stats.terms(:,dim), 2) == sum (stats.terms, 2));
  Nb = 1 + sum(df(k));
  Nt = numel (k);
  L = zeros (n, sum (df) + 1);
  for j = 1:Nt
    L(:, i(k(j)) - df(k(j)) + 1 : i(k(j))) = stats.X(:,i(k(j)) - ...
                                             df(k(j)) + 1 : i(k(j)));
  end
  L(:,1) = 1;
  H = unique (L, 'rows', 'stable');
  Ng = size (H, 1);
  idx = zeros (Ng, 1);
  for k = 1:Ng
    idx(k) = find (all (L == H(k, :), 2),1);
  end

  % Define bootfun for bootstraping the model residuals and returning the group means
  bootfun = @(r) H * lmfit (X, fitted + r .* se, W);

  % Perform bootstrap
  if nargout > 0
    warning ('off','bootknife:lastwarn')
    coeffs = bootknife (resid, nboot, bootfun, alpha, [], ncpus);
    warning ('on','bootknife:lastwarn')
  else
    coeffs = [];
    bootknife (resid, nboot, bootfun, alpha, [], ncpus);
  end

end
