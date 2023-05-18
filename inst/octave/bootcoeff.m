% -- Function File: bootcoeff (STATS)
% -- Function File: bootcoeff (STATS, NBOOT)
% -- Function File: bootcoeff (STATS, NBOOT, PROB)
% -- Function File: bootcoeff (STATS, NBOOT, PROB)
% -- Function File: bootcoeff (STATS, NBOOT, PROB, SEED)
% -- Function File: bootcoeff (STATS, NBOOT, PROB, PRIOR, SEED)
% -- Function File: COEFF = bootcoeff (STATS, ...)
% -- Function File: [COEFF, BOOTSTAT] = bootcoeff (STATS, ...)
%
%     'bootcoeff (STATS)' uses the STATS structure output from fitlm or anovan
%     functions (from the v1.5+ of the Statistics package in OCTAVE) and
%     Bayesian nonparametric bootstrap [1] to compute and return the following
%     statistics:
%        • original: the coefficient(s) from regressing y on X
%        • bias: bootstrap bias estimate(s)
%        • median: the median of the posterior distribution(s)
%        • CI_lower: lower bound(s) of the 95% credible interval
%        • CI_upper: upper bound(s) of the 95% credible interval
%        • p-val: two-tailed p-value(s) for the parameter(s) being equal to 0
%          By default, the credible intervals are shortest probability intervals,
%          which represent a more computationally stable version of the highest
%          posterior density interval [2]. The p-value(s) is/are computed from
%          the Student-t (null) distribution(s) constructed from the posterior
%          statistics and heteroscedasticity-consistent standard errors [3,4].
%
%     'bootcoeff (STATS, NBOOT)' specifies the number of bootstrap resamples,
%     where NBOOT must be a positive integer. If empty, the default value of
%     NBOOT is the scalar: 2000.
%
%     'bootcoeff (STATS, NBOOT, PROB)' where PROB is numeric and sets the lower
%     lower and upper bounds of the credible interval(s). The value(s) of
%     PROB must be between 0 and 1. PROB can either be:
%        • scalar: To set the central mass of shortest probability intervals
%                  (SPI) to 100*(1-PROB)%
%        • vector: A pair of probabilities defining the lower and upper
%                  percentiles of the credible interval(s) as 100*(PROB(1))%
%                  and 100*(PROB(2))% respectively. 
%          Credible intervals are not calculated when the value(s) of PROB
%          is/are NaN. The default value of PROB is the scalar 0.95.
%
%     'bootcoeff (STATS, NBOOT, PROB, PRIOR)' accepts a positive real numeric
%     scalar to parametrize the form of the symmetric Dirichlet distribution.
%     The Dirichlet distribution is the conjugate PRIOR used to randomly
%     generate weights for linear least squares fitting to the observed data,
%     and subsequently to estimate the posterior for the regression coefficients
%     by Bayesian bootstrap. If PRIOR is not provided, or is empty, it will be
%     set to 1, corresponding to Bayes rule: a uniform (or flat) Dirichlet
%     distribution (in the range [0, 1]). For a stronger prior, set PRIOR to 
%     > 1, for example, if the experiment is confirmatory and reproduced. For a
%     weaker prior, set PRIOR to < 1 (e.g. 0.5 for Jeffrey's prior). Jeffrey's
%     prior may be appropriate for estimates from small samples (n < 10) in 
%     exploratory settings, where the amount of data may be assumed to
%     inadequately define the parameter space. 
%
%     'bootcoeff (STATS, NBOOT, PROB, PRIOR, SEED)' initialises the Mersenne
%     Twister random number generator using an integer SEED value so that
%     'bootcoeff' results are reproducible.
%
%     'COEFF = bootcoeff (STATS, ...) returns a structure with the following
%     fields (defined above): original, bias, median, CI_lower, CI_upper, tstat
%     and pval. These statistics summarise the posterior distributions of the
%     coefficients from the linear model.
%
%     '[COEFF, BOOTSTAT] = bootcoeff (STATS, ...) also returns the bootstrap
%     statistics (i.e. posterior) for the coefficients.
%
%  Requirements: bootcoeff is only supported in GNU Octave and requires the
%  Statistics package version 1.5+.
%
%  See also `bootbayes`.
%
%  Bibliography:
%  [1] Rubin (1981) The Bayesian Bootstrap. Ann. Statist. 9(1):130-134
%  [2] Liu, Gelman & Zheng (2015). Simulation-efficient shortest probability
%        intervals. Statistics and Computing, 25(4), 809–819.
%  [3] Hall and Wilson (1991) Two Guidelines for Bootstrap Hypothesis Testing.
%        Biometrics, 47(2), 757-762
%  [4] Long and Ervin (2000) Using Heteroscedasticity Consistent Standard
%        Errors in the Linear Regression Model. Am. Stat, 54(3), 217-224
%
%  bootcoeff (version 2023.05.18)
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


function [coeffs, bootstat] = bootcoeff (STATS, nboot, prob, prior, seed)

  % Check input aruments
  if (nargin < 1)
    error ('bootcoeff usage: ''bootcoeff (STATS)'' atleast 1 input arguments required');
  end
  if (nargin < 2)
    nboot = 2000;
  end
  if (nargin < 3)
    prob = 0.95;
  end
  if (nargin < 4)
    prior = []; % Use default in bootbayes
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
  if (nargin > 4)
    if (ISOCTAVE)
      randg ('seed', seed);
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
      bootbayes (y, X, nboot, prob, prior);
    case 1
      coeffs = bootbayes (y, X, nboot, prob, prior);
    otherwise
      [coeffs, bootstat] = bootbayes (y, X, nboot, prob, prior);
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
%! [P,ATAB,STATS] = anovan (dv,g,'contrasts','treatment');
%! STATS.coeffnames
%! # Uniform prior, 95% credible intervals
%! bootcoeff (STATS,2000,0.95,1.0)