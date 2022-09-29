%  Function File: bootknife
%
%  coeffs = bootcoeff(STATS)
%
%  Bootstrap the regression coefficients from a linear model. bootcoeff accepts
%  as input the STATS structure from fitlm or anovan and returns a structure,
%  coeffs, which contains the following fields:
%    original: contains the regression coefficients from the original data
%    bias: contains the bootstrap estimate of bias
%    std_error: contains the bootstrap standard error
%    CI_lower: contains the lower bound of the bootstrap confidence interval
%    CI_upper: contains the upper bound of the bootstrap confidence interval
%  The confidence intervals are 95% bias-corrected intervals. If no output is
%  requested, the results are printed to stdout.
%
%  bootcoef is only supported in GNU Octave and requires the Statistics package
%  version 1.5 or later.
%
%  bootcoeff v0.1.0.0 (27/09/2022)
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


function coeffs = bootcoeff (STATS, nboot, alpha, ncpus)

  % Error checking
  if (nargin < 2)
    nboot = 2000;
  end
  if (nargin < 3)
    alpha = 0.05;
  end
  if (nargin < 2)
    ncpus = 0;
  end
  % Check if running in Octave (else assume Matlab)
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));
  if ~ISOCTAVE
    error ('bootcoeff is only supported in by Octave')
  end
  statspackage = ismember ({info.Name}, 'statistics');
  if any (statspackage)
    if str2num(info(statspackage).Version(1:3)) < 1.5
      error ('bootcoeff requires version >= 1.5 of the statistics package')
    end
  end
    
  % Fetch required information from STATS structure
  X = STATS.X;
  b = STATS.coeffs(:,1);
  fitted = X * b;
  lmfit = STATS.lmfit;
  W = full (STATS.W);
  se = sqrt (diag (W));
  resid = STATS.resid;

  % Define bootfun for bootstraping the model residuals
  bootfun = @(r) lmfit (X, fitted + r./se, W);

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