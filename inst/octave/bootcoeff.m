% -- Function File: bootcoeff (STATS)
% -- Function File: bootcoeff (STATS, CLUSTID)
% -- Function File: bootcoeff (STATS, BLOCKSZ)
% -- Function File: bootcoeff (STATS, ..., NBOOT)
% -- Function File: bootcoeff (STATS, ..., NBOOT, PROB)
% -- Function File: bootcoeff (STATS, ..., NBOOT, PROB)
% -- Function File: bootcoeff (STATS, ..., NBOOT, PROB, SEED)
% -- Function File: bootcoeff (STATS, ..., NBOOT, PROB, PRIOR, SEED)
% -- Function File: COEFF = bootcoeff (STATS, ...)
% -- Function File: [COEFF, BOOTSTAT] = bootcoeff (STATS, ...)
%
%     'bootcoeff (STATS)' accepts the STATS structure output from fitlm or anovan
%     functions (from the v1.5+ of the Statistics package in OCTAVE) and uses
%     the `bootbayes` and `bootwild` functions to compute and print 95% credible
%     intervals and p-values respectively for each model coefficient. Please see
%     the help for both of these functions for more information.
%
%     'bootcoeff (STATS, CLUSTID)' specifies a vector or cell array of numbers
%     or strings respectively to be used as cluster labels or identifiers.
%     Rows of the data with the same CLUSTID value are treated as clusters with
%     dependent errors. If empty (default), no clustered resampling is performed
%     and all errors are treated as independent.
%
%     'bootcoeff (STATS, BLOCKSZ)' specifies a scalar, which sets the block size
%     for bootstrapping when the errors have serial dependence. Rows of the data
%     within the same block are treated as having dependent errors. If empty
%     (default), no clustered resampling is performed and all errors are treated
%     as independent.
%
%     'bootcoeff (STATS, ..., NBOOT)' specifies the number of bootstrap
%     resamples, where NBOOT must be a positive integer. If empty, the default
%     value of NBOOT is the scalar: 2000.
%
%     'bootcoeff (STATS, ..., NBOOT, PROB)' where PROB is numeric and sets
%     the lower and upper bounds of the credible interval(s). The value(s) of
%     PROB must be between 0 and 1. PROB can either be:
%        • scalar: To set the central mass of shortest probability intervals
%                  (SPI) to 100*(1-PROB)%
%        • vector: A pair of probabilities defining the lower and upper
%                  percentiles of the credible interval(s) as 100*(PROB(1))%
%                  and 100*(PROB(2))% respectively. 
%          Credible intervals are not calculated when the value(s) of PROB
%          is/are NaN. The default value of PROB is the scalar 0.95.
%
%     'bootcoeff (STATS, ..., NBOOT, PROB, PRIOR)' accepts a positive real
%     numeric scalar to parametrize the form of the symmetric Dirichlet
%     distribution. The Dirichlet distribution is the conjugate PRIOR used to
%     randomly generate weights for linear least squares fitting of the observed
%     data, and subsequently to estimate the posterior for the regression
%     coefficients by Bayesian bootstrap. If PRIOR is not provided, or is empty,
%     it will be set to 1, corresponding to Bayes rule: a uniform (or flat)
%     Dirichlet distribution (in the range [0, 1]). For a weaker prior, set
%     PRIOR to < 1 (e.g. 0.5 for Jeffrey's prior).
%
%     'bootcoeff (STATS, ..., NBOOT, PROB, PRIOR, SEED)' initialises the
%     Mersenne Twister random number generator using an integer SEED value so
%     that `bootcoeff` results are reproducible.
%
%     'COEFF = bootcoeff (STATS, ...) returns a structure with the following
%     fields (defined above): original, bias, median, CI_lower, CI_upper, tstat,
%     pval and fpr.
%
%     '[COEFF, BOOTSTAT] = bootcoeff (STATS, ...) also returns the bootstrap
%     statistics for the coefficients.
%
%  Requirements: bootcoeff is only supported in GNU Octave and requires the
%  Statistics package version 1.5+.
%
%  See also `bootbayes` and `bootwild`.
%
%  bootcoeff (version 2023.06.07)
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


function [coeffs, bootstat] = bootcoeff (STATS, arg2, nboot, prob, prior, seed)

  % Check input aruments
  if (nargin < 1)
    error ('bootcoeff usage: ''bootcoeff (STATS)'' atleast 1 input arguments required');
  end
  if (nargin < 2)
    arg2 = [];
  end
  if (nargin < 3)
    nboot = 2000;
  end
  if (nargin < 4)
    prob = 0.95;
  end
  if (nargin < 5)
    prior = []; % Use default in bootbayes
  end
  if (nargin < 6)
    seed = []; % Use default in bootbayes
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
  if (nargin > 5)
    if (ISOCTAVE)
      randg ('seed', seed);
      randn ('seed', seed);
    else
      rng ('default');
    end
  end

  % Fetch required information from STATS structure
  X = full (STATS.X);
  b = STATS.coeffs(:,1);
  fitted = X * b;
  resid = STATS.resid;
  y = fitted + resid;
  N = numel (resid);
  if (~ all (diag (full (STATS.W) == 1)))
    error ('bootcoeff: Incompatible with the ''weights'' argument in ''anovan'' or ''fitlm''')
  end

  % Perform Bayesian bootstrap
  switch (nargout)
    case 0
      bootbayes (y, X, arg2, nboot, prob, prior);
      bootwild (y, X, arg2, nboot);
    case 1
      coeffs = bootbayes (y, X, arg2, nboot, prob, prior);
      nhst = bootwild (y, X, arg2, nboot);
      coeffs.tstat = nhst.tstat;
      coeffs.pval = nhst.pval;
    otherwise
      [coeffs, bootstat] = bootbayes (y, X, arg2, nboot, prob, prior);
      nhst = bootwild (y, X, arg2, nboot);
      coeffs.tstat = nhst.tstat;
      coeffs.pval = nhst.pval;
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
%! [P,ATAB,STATS] = anovan (dv,g,'contrasts','treatment','sstype',2);
%! STATS.coeffnames
%! # Uniform prior, 95% credible intervals
%! bootcoeff (STATS,[],2000,0.95,1.0)