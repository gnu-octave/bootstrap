% -- Function File: bootemm (STATS, DIM)
% -- Function File: bootemm (STATS, DIM, CLUSTID)
% -- Function File: bootemm (STATS, DIM, BLOCKSZ)
% -- Function File: bootemm (STATS, DIM, ..., NBOOT)
% -- Function File: bootemm (STATS, DIM, ..., NBOOT, PROB)
% -- Function File: bootemm (STATS, DIM, ..., NBOOT, PROB)
% -- Function File: bootemm (STATS, DIM, ..., NBOOT, PROB, SEED)
% -- Function File: EMM = bootemm (STATS, DIM, ...)
% -- Function File: [EMM, BOOTSTAT] = bootemm (STATS, DIM, ...)
%
%     'bootemm (STATS, DIM)' accepts the STATS structure output from fitlm or
%     anovan functions (from the v1.5+ of the Statistics package in OCTAVE)
%     and uses the `bootbayes` function to compute and print 95% credible
%     intervals each estimated marginal mean along the dimension DIM. Please
%     see the help documentation in `bootbayes` for more information.
%
%     'bootemm (STATS, DIM, CLUSTID)' specifies a vector or cell array of
%     numbers or strings respectively to be used as cluster labels or
%     identifiers. Rows of the data with the same CLUSTID value are treated
%     as clusters with dependent errors. If empty (default), no clustered
%     resampling is performed and all errors are treated as independent.
%
%     'bootemm (STATS, DIM, BLOCKSZ)' specifies a scalar, which sets the block
%     size for bootstrapping when the errors have serial dependence. Rows of
%     the data within the same block are treated as having dependent errors. If
%     empty (default), no clustered resampling is performed and all errors are
%     treated as independent.
%
%     'bootemm (STATS, DIM, ..., NBOOT)' specifies the number of bootstrap
%     resamples, where NBOOT must be a positive integer. If empty, tHe default
%     value of NBOOT is the scalar: 2000.
%
%     'bootemm (STATS, DIM, ..., NBOOT, PROB)' where PROB is numeric and
%     sets the lower and upper bounds of the credible interval(s). The value(s)
%     of PROB must be between 0 and 1. PROB can either be:
%        • scalar: To set the central mass of shortest probability intervals
%                  (SPI) to 100*(1-PROB)%
%        • vector: A pair of probabilities defining the lower and upper
%                  percentiles of the credible interval(s) as 100*(PROB(1))%
%                  and 100*(PROB(2))% respectively. 
%          Credible intervals are not calculated when the value(s) of PROB
%          is/are NaN. The default value of PROB is the scalar 0.95.
%
%     'bootemm (STATS, DIM, ..., NBOOT, PROB, PRIOR)' accepts a positive
%     real numeric scalar to parametrize the form of the symmetric Dirichlet
%     distribution. The Dirichlet distribution is the conjugate PRIOR used to
%     randomly generate weights for linear least squares fitting of the observed
%     data, and subsequently to estimate the posterior for the regression
%     coefficients by Bayesian bootstrap. If PRIOR is not provided, or is empty,
%     it will be set to 1, corresponding to Bayes rule: a uniform (or flat)
%     Dirichlet distribution (in the range [0, 1]). For a weaker prior, set
%     PRIOR to < 1 (e.g. 0.5 for Jeffrey's prior).
%
%     'bootemm (STATS, DIM, ..., NBOOT, PROB, PRIOR, SEED)', initialises the
%     Mersenne Twister random number generator using an integer SEED value so
%     that 'bootemm' results are reproducible.
%
%     'EMM = bootemm (STATS, DIM, ...) returns a structure with the following
%     fields (defined above): original, bias, median, CI_lower & CI_upper.
%
%     '[EMM, BOOTSTAT] = bootemm (STATS, DIM, ...) also returns the bootstrap
%     statistics for the estimated marginal means.
%
%  See also `bootbayes`.
%
%  Requirements: bootemm is only supported in GNU Octave and requires the
%  Statistics package version 1.5+.
%
%  bootemm (version 2023.06.07)
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


function [emm, bootstat] = bootemm (STATS, dim, dep, nboot, prob, prior, seed)

  % Check input aruments
  if (nargin < 2)
    error ('bootemm usage: ''bootemm (STATS, dim)'' atleast 2 input arguments required');
  end
  if (nargin < 3)
    dep = [];
  end
  if (nargin < 4)
    nboot = 2000;
  end
  if (nargin < 5)
    prob = 0.95;
  end
  if (nargin < 6)
    prior = [];
  end
  if (nargin < 7)
    seed = []; % Use default in bootbayes
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
  N = numel (STATS.contrasts);
  for j = 1:N
    if (isnumeric (STATS.contrasts{j}))
      % Check that the columns sum to 0
      if (any (abs (sum (STATS.contrasts{j})) > eps("single")))
         error (strcat(["Use a STATS structure from a model refit with"], ...
                       [" sum-to-zero contrast coding, e.g. ""simple"""]));
      end
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

  % Prepare the hypothesis matrix (H)
  df = STATS.df;
  i = 1 + cumsum(df);
  k = find (sum (STATS.terms(:,dim), 2) == sum (STATS.terms, 2));
  Nt = numel (k);
  H = zeros (N, sum (df) + 1);
  for j = 1:Nt
    H(:, i(k(j)) - df(k(j)) + 1 : i(k(j))) = STATS.X(:,i(k(j)) - ...
                                             df(k(j)) + 1 : i(k(j)));
  end
  H(:,1) = 1;
  L = unique (H, 'rows', 'stable');
  Ng = size (L, 1);
  idx = zeros (Ng, 1);
  for k = 1:Ng
    idx(k) = find (all (H == L(k, :), 2),1);
  end

  % Perform Bayesian bootstrap
  switch (nargout)
    case 0
      bootbayes (y, X, dep, nboot, prob, prior, seed, L);
    case 1
      emm = bootbayes (y, X, dep, nboot, prob, prior, seed, L);
    otherwise
      [emm, bootstat] = bootbayes (y, X, dep, nboot, prob, prior, seed, L);
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
%! # Uniform prior, 95% credible intervals
%! bootemm (STATS, DIM, [], 2000, 0.95, 1.0); 
