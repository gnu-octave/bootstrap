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
%  Semi-parametric bootstrap of the estimated marginal means from a linear model.
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
%  [2,3]. By default, the confidence intervals constructed are expanded bias-
%  corrected percentile intervals [4,5]. The list of means and their bootstrap
%  statistics correspond to the names STATS.grpnames. If no output is requested,
%  the results are printed to stdout*. 
%
%  * Note that the confidence interval coverage printed reflects the coverage
%  after expansion using Student's t-distribution (if applicable) [5].When
%  reporting confidence intervals, the nominal coverage should be used, 
%  e.g. 95% if alpha is 0.05.
%
%  EMMEANS = bootemm (STATS, NBOOT) also specifies the number of bootstrap
%  samples. NBOOT must be a scalar. By default, NBOOT is 2000.
%
%  EMMEANS = bootemm (STATS, NBOOT, ALPHA) where ALPHA is numeric and
%  sets the lower and upper bounds of the confidence interval(s). The value(s)
%  of ALPHA must be between 0 and 1. ALPHA can either be:
%
%  1) a scalar value to set the (nominal) central coverage to 100*(1-ALPHA)%
%  with (nominal) lower and upper percentiles of the confidence intervals at
%  100*(ALPHA/2)% and 100*(1-ALPHA/2)% respectively. The intervals constructed
%  are expanded [5] and bias-corrected [4] to improve central coverage and
%  accuracy.
%
%  2) a vector containing a pair of quantiles to set the (nominal) lower and
%  upper percentiles of the confidence interval(s) at 100*(ALPHA(1))% and
%  100*(ALPHA(2))%. The intervals constructed are simple percentile intervals.
%  
%  Confidence interval endpoints are not calculated when the value(s) of ALPHA
%  is/are NaN. If empty (or not specified), the default value for ALPHA is 0.05
%
%  EMMEANS = bootemm (STATS, NBOOT, ALPHA, NPROC) also sets the number of
%  parallel processes to use to accelerate computations on multicore machines.
%  This feature requires the Parallel package (in Octave).
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
%  [4] Efron (1981) Nonparametric Standard Errors and Confidence Intervals.
%        Can J Stat. 9(2:139-158
%  [5] Hesterberg, Tim (2014), What Teachers Should Know about the 
%        Bootstrap: Resampling in the Undergraduate Statistics Curriculum, 
%        http://arxiv.org/abs/1411.5279
%
%  bootemm (version 2022.12.16)
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
    boot (1, 1, false, seed);
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
  dfe = stats.dfe;

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

  % Perform adjustments to improve coverage for small samples
  switch (numel (alpha))
    case 1
      % Expanded bias-corrected confidence intervals (recommended)
      % Create distribution functions
      stdnormcdf = @(x) 0.5 * (1 + erf (x / sqrt (2)));
      % If bootfun is the mean, adjust alpha level to expand percentiles
      % using Student's t-distribution for improved central coverage in
      % when sample size is small
      if exist('betaincinv','file')
        studinv = @(p, df) - sqrt ( df ./ betaincinv (2 * p, df / 2, 0.5) - df);
      else
        % Earlier versions of matlab do not have betaincinv
        % Instead, use betainv from the Statistics and Machine Learning Toolbox
        try 
          studinv = @(p, df) - sqrt ( df ./ betainv (2 * p, df / 2, 0.5) - df);
        catch
          % Use the Normal distribution (i.e. do not expand quantiles) if
          % either betaincinv or betainv are not available
          studinv = @(p,df) sqrt (2) * erfinv (2 * p-1);
        end
      end
      adj_alpha = stdnormcdf (studinv (alpha / 2, stats.dfe)) * 2;
    case 2
      % Simple percentile confidence intervals (no adjustment)
      % No adjustment
      adj_alpha = alpha;
  end
  % Check that the tail probabilities are > 0 and < 1
  if ~ all (adj_alpha > 0 & adj_alpha < 1)
    error ('bootknife: the sample size is too small for the specified ALPHA')
  end

  % Perform bootstrap
  if nargout > 0
    warning ('off','bootknife:lastwarn')
    emmeans = bootknife (resid, nboot, bootfun, adj_alpha, [], ncpus);
    warning ('on','bootknife:lastwarn')
  else
    bootknife (resid, nboot, bootfun, adj_alpha, [], ncpus);
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