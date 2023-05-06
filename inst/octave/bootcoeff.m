% -- Function File: bootcoeff (STATS)
% -- Function File: bootcoeff (STATS, NBOOT)
% -- Function File: bootcoeff (STATS, NBOOT, ALPHA)
% -- Function File: bootcoeff (STATS, NBOOT, ALPHA)
% -- Function File: bootcoeff (STATS, NBOOT, ALPHA, SEED)
% -- Function File: COEFF = bootcoeff (STATS, ...)
% -- Function File: [COEFF, BOOTSTAT] = bootcoeff (STATS, ...)
%
%     'bootcoeff (STATS)' uses the STATS structure output from fitlm or anovan
%     functions (from the v1.5+ of the Statistics package in OCTAVE) and
%     Bayesian bootstrap to compute and return the following statistics:
%        • original: regression coefficients from the original data
%        • bias: bootstrap estimate of the bias of the coefficients
%        • std_error: bootstrap estimate of the standard error
%        • CI_lower: lower bound of the 95% bootstrap confidence interval
%        • CI_upper: upper bound of the 95% bootstrap confidence interval
%          Here, the confidence intervals, or credible intervals in the context
%          of the Bayesian statistical framework, are percentile intervals [2].
%
%     'bootcoeff (STATS, NBOOT)' specifies the number of bootstrap resamples,
%     where NBOOT must be a positive integer. If empty, tHe default value of
%     NBOOT is the scalar: 2000.
%
%     'bootcoeff (STATS, NBOOT, ALPHA)', where ALPHA is numeric and sets the
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
%     'bootcoeff (STATS, NBOOT, ALPHA, SEED)' initialises the Mersenne
%     Twister random number generator using an integer SEED value so that
%     bootcoeff results are reproducible.
%
%     'COEFF = bootcoeff (STATS, ...) returns a structure with the following
%     fields (defined above): original, bias, std_error, CI_lower, CI_upper.
%     These statistics correspond to the coefficients of the linear model.
%
%     '[COEFF, BOOTSTAT] = bootcoeff (STATS, ...) also returns the bootstrap
%     statistics for the coefficients.
%
%  Requirements: bootcoeff is only supported in GNU Octave and requires the
%  Statistics package version 1.5+.
%
%  Bibliography:
%  [1] Rubin (1981) The Bayesian Bootstrap. Ann. Statist. 9(1):130-134
%  [2] Efron and Tibshirani. Chapter 16 Hypothesis testing with the
%       bootstrap in An introduction to the bootstrap (CRC Press, 1994)
%
%  bootcoeff (version 2023.05.02)
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


function [coeffs, bootstat] = bootcoeff (STATS, nboot, alpha, seed)

  % Check input aruments
  if (nargin < 1)
    error ('bootcoeff usage: ''bootcoeff (STATS)'' atleast 1 input arguments required');
  end
  if (nargin < 2)
    nboot = 2000;
  end
  if (nargin < 3)
    alpha = [0.025, 0.975];
  end
  if (nargin > 3)
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
    error ('bootcoeff: Only supported by Octave')
  end
  statspackage = ismember ({info.Name}, 'statistics');
  if ((~ any (statspackage)) || (str2double (info (statspackage).Version(1:3)) < 1.5))
    error ('bootcoeff: Requires version >= 1.5 of the statistics package')
  end

  % Fetch required information from STATS structure
  X = full (STATS.X);
  b = STATS.coeffs(:,1);
  fitted = X * b;
  resid = STATS.resid;
  y = fitted + resid;
  if (~ all (diag (full (STATS.W) == 1)))
    error ('bootcoeff: Incompatible with the ''weights'' argument in ''anovan'' or ''fitlm''')
  end

  % Perform Bayesian bootstrap
  switch (nargout)
    case 0
      bootbayes (y, X, nboot, alpha);
    case 1
      coeffs = bootbayes (y, X, nboot, alpha);
    otherwise
      [coeffs, bootstat] = bootbayes (y, X, nboot, alpha);
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