% -- Function File: bootwild (y)
% -- Function File: bootwild (y, X)
% -- Function File: bootwild (y, X, CLUSTID)
% -- Function File: bootwild (y, X, BLOCKSZ)
% -- Function File: bootwild (y, X, .., NBOOT)
% -- Function File: bootwild (y, X, .., NBOOT, SEED)
% -- Function File: STATS = bootwild (y, ...)
% -- Function File: [STATS, BOOTSTAT] = bootwild (y, ...)
%
%     'bootwild (y)' performs a null hypothesis significance test for the
%     mean of y being equal to 0. This function implements wild bootstrap-t
%     resampling of the Rademacher distribution of the residuals, and computes
%     p-values after imposing the null hypothesis (H0) [1,2]. The following
%     statistics are printed to the standard output:
%        • original: the mean of the data vector y
%        • tstat: bootstrap bias estimate(s)
%        • pval: two-tailed p-value(s) for the parameter(s) being equal to 0
%        • fpr: minimum false positive risk for the corresponding p-value
%          The p-values are computed following both of the guidelines by Hall
%          and Wilson [3]. The minimum false positive risk (FPR) is computed
%          according to the Sellke-Berger approach as described in [4,5].
%
%     'bootwild (y, X)' also specifies the design matrix (X) for least squares
%     regression of y on X. X should be a column vector or matrix the same
%     number of rows as y. If the X input argument is empty, the default for X
%     is a column of ones (i.e. intercept only) and thus the statistic computed
%     reduces to the mean (as above). The statistics calculated and returned in
%     the output then relate to the coefficients from the regression of y on X.
%
%     'bootwild (y, X, CLUSTID)' specifies a vector or cell array of numbers
%     or strings respectively to be used as cluster labels or identifiers.
%     Rows in y (and X) with the same CLUSTID value are treated as clusters
%     with dependent errors. Rows of y (and X) assigned to a particular
%     cluster will have identical sign-flipping during wild bootstrap. If empty
%     (default), no clustered resampling is performed and all errors are
%     treated as independent.
%
%     'bootwild (y, X, BLOCKSZ)' specifies a scalar, which sets the block size
%     for bootstrapping when the residuals have serial dependence. Identical
%     sign-flipping occurs within each (consecutive) block of length BLOCKSZ
%     during wild bootstrap. Rows of y (and X) within the same block are
%     treated as having dependent errors. If empty (default), no block
%     resampling is performed and all errors are treated as independent.
%
%     'bootwild (y, X, ..., NBOOT)' specifies the number of bootstrap resamples,
%     where NBOOT must be a positive integer. If empty, the default value of
%     NBOOT is 2000.
%
%     'bootwild (y, X, ..., NBOOT, SEED)' initialises the Mersenne Twister
%     random number generator using an integer SEED value so that 'bootwild'
%     results are reproducible.
%
%     'STATS = bootwild (...) returns a structure with the following fields
%     (defined above): original, tstat, pval and fpr. 
%
%     '[STATS, BOOTSTAT] = bootwild (...)  also returns the a vector (or
%     matrix) of bootstrap statistics (BOOTSTAT) calculated over the bootstrap
%     resamples.
%
%  Bibliography:
%  [1] Wu (1986). Jackknife, bootstrap and other resampling methods in
%        regression analysis (with discussions). Ann Stat.. 14: 1261–1350. 
%  [2] Cameron, Gelbach and Miller (2008) Bootstrap-based Improvements for
%        Inference with Clustered Errors. Rev Econ Stat. 90(3), 414-427
%  [3] Hall and Wilson (1991) Two Guidelines for Bootstrap Hypothesis Testing.
%        Biometrics, 47(2), 757-762
%  [4] Colquhoun (2019) The False Positive Risk: A Proposal Concerning What
%        to Do About p-Values, Am Stat. 73:sup1, 192-201
%  [5] Sellke, Bayarri and Berger (2001) Calibration of p-values for Testing
%        Precise Null Hypotheses. Am Stat. 55(1), 62-71
%
%  bootwild (version 2023.06.07)
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


function [stats, bootstat] = bootwild (y, X, arg3, nboot, seed)

  % Check the number of function arguments
  if (nargin < 1)
    error ('bootwild: y must be provided');
  end
  if (nargin > 8)
    error ('bootwild: Too many input arguments')
  end
  if (nargout > 2)
    error ('bootwild: Too many output arguments')
  end

  % Check if running in Octave (else assume Matlab)
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

  % Calculate the length of y
  if (nargin < 1)
    error ('bootwild: DATA must be provided');
  end
  sz = size (y);
  if ( (sz(1) < 2) || (sz (2) > 1) )
    error ('bootwild: y must be a column vector');
  end
  n = numel (y);

  % Evaluate the design matrix
  if ( (nargin < 2) || (isempty (X)) )
    X = ones (n, 1);
  end

  % Calculate number of parameters
  p = size (X, 2);

  % Evaluate cluster IDs or block size
  if ( (nargin > 2) && (~ isempty (arg3)) )
    if (isscalar (arg3))
      % Prepare for wild block bootstrap
      blocksz = arg3;
      G = fix (n / blocksz);
      IC = (G + 1) * ones (n, 1);
      IC(1 : blocksz * G, :) = reshape (ones (blocksz, 1) * [1 : G], [], 1);
      G = IC(end);
      method = 'block ';
    else
      % Prepare for wild cluster bootstrap
      clustid = arg3;
      if (bsxfun (@ne, size (clustid), sz))
        error ('bootwild: clustid must be the same size as y')
      end
      [C, IA, IC] = unique (clustid);
      G = numel (C); % Number of clusters
      method = 'cluster ';
    end
    UC = unique (IC);
    clusters = cell2mat (cellfun (@(i) IC == UC(i), ...
                                  num2cell (1:G), 'UniformOutput', false));
  else
    G = n;
    IC = [];
    clusters = [];
    method = "";
  end

  % Evaluate number of bootstrap resamples
  if ( (nargin < 4) || (isempty (nboot)) )
    nboot = 2000;
  else
    if (~ isa (nboot, 'numeric'))
      error ('bootwild: NBOOT must be numeric');
    end
    if (numel (nboot) > 1)
      error ('bootwild: NBOOT must be scalar');
    end
    if (nboot ~= abs (fix (nboot)))
      error ('bootwild: NBOOT must be a positive integers');
    end
  end

  % Set random seed
  if ( (nargin > 4) && (~ isempty (seed)) )
    if (ISOCTAVE)
      randn ('seed', seed);
    else
      rng (seed);
    end
  end

  % Create least squares anonymous function for bootstrap
  bootfun = @(y) lmfit (X, y, clusters);

  % Calculate estimate(s)
  S = bootfun (y);
  original = S.b;
  std_err = S.se;
  t = original ./ std_err;

  % Perform sign flipping of the residuals and create resamples (Rademacher distribution)
  s = sign (randn (G, nboot));
  if (~ isempty (IC))
    s = s(IC, :);  % Enforce clustering/blocking
  end
  yf = X * original;
  r = y - yf;
  Y = bsxfun (@plus, yf, r .* s);

  % Compute bootstap statistics
  bootout = cell2mat (cellfun (bootfun, num2cell (Y, 1), 'UniformOutput', false));
  bootstat = [bootout.b];
  bootse = [bootout.se];

  % Enforce H0 and studentize the bootstrap statistics following both
  % guidelines described in Hall and Wilson (1991) Biometrics, 47(2), 757-762
  T = bsxfun (@minus, bootstat, original) ./ bootse;

  % Compute two-tailed p-values
  pval = nan (p, 1);
  for j = 1:p
    if ( isnan (std_err(j)) )
      pval(j) = NaN;
    else
      [cdf, x] = empcdf (abs (T(j,:)));
      pval(j) = 1 - interp1 (x, cdf, abs (t(j)), 'linear', 1);
    end
  end

  % Compute minimum false positive risk
  fpr = pval2fpr (pval);

  % Prepare output arguments
  stats = struct;
  stats.original = original;
  stats.std_err = std_err;
  stats.tstat = t;
  stats.pval = pval;
  stats.fpr = fpr;

  % Print output if no output arguments are requested
  if (nargout == 0) 
    print_output (stats, nboot, p, method);
  end

end

%--------------------------------------------------------------------------

%% FUNCTION TO FIT THE LINEAR MODEL

function S = lmfit (X, y, clusters)

  % Get model coefficients by solving the linear equation by matrix arithmetic

  % Solve linear equation to minimize least squares and compute the
  % regression coefficients (b)
  invG = pinv (X' * X);
  b = invG * (X' * y);

  % Calculate heteroscedasticity-consistent (HC) or cluster robust (CR) standard 
  % errors (CR) for the regression coefficients. When the number of observations
  % equals the number of clusters, the calculations for CR reduce to HC.
  % References: 
  %   Long and Ervin (2000) Am. Stat, 54(3), 217-224
  %   Cameron, Gelbach and Miller (2008) Rev Econ Stat. 90(3), 414-427
  %   MacKinnon & Webb (2020) QED Working Paper Number 1421
  yf = X * b;
  u = y - yf;
  if ( (nargin < 3) || isempty (clusters) )
    % Heteroscedasticity-Consistent (HC0) standard errors
    meat = X' * diag (u.^2) * X;
  else
    % Cluster Robust (CR0) standard errors
    Sigma = cellfun (@(g) X(g,:)' * u(g) * u(g)' * X(g,:), ...
                  num2cell (clusters, 1), 'UniformOutput', false);
    meat = sum (cat (3, Sigma{:}), 3);
  end
  se = sqrt (diag (invG * meat * invG));

  % Prepare output
  S.b = b;
  S.se = se;

end

%--------------------------------------------------------------------------

%% FUNCTION TO COMPUTE EMPIRICAL DISTRIBUTION FUNCTION

function [F, x] = empcdf (y)

  % Subfunction to calculate empirical cumulative distribution function

  % Check input argument
  if (~ isa (y, 'numeric'))
    error ('bootwild:empcdf: y must be numeric');
  end
  if (all (size (y) > 1))
    error ('bootwild:empcdf: y must be a vector');
  end
  if (size (y, 2) > 1)
    y = y.';
  end

  % Discard NaN values
  ridx = isnan (y);
  y(ridx) = [];

  % Get size of y
  N = numel (y);

  % Create empirical CDF
  x = sort (y);
  [x, F] = unique (x,'rows','last');
  F = (F - 1) / (N - 1);

end

%--------------------------------------------------------------------------

%% FUNCTION TO COMPUTE FALSE POSITIVE RISK (FPR)

function fpr = pval2fpr (p)

  % Subfunction to compute minimum false positive risk. These are calculated
  % from a Bayes factor based on the sampling distributions of the p-value and
  % that H0 and H1 have equal prior probabilities. This is called the Sellke-
  % Berger approach.
  % 
  % References:
  %  Held and Ott (2018) On p-Values and Bayes Factors. 
  %    Annu. Rev. of Stat. Appl. 5:393-419
  %  David Colquhoun (2019) The False Positive Risk: A Proposal 
  %    Concerning What to Do About p-Values, The American Statistician, 
  %    73:sup1, 192-201, DOI: 10.1080/00031305.2018.1529622 

  % Calculate minimum Bayes Factor (P(H0) / P(H1)) by the Sellke-Berger method 
  logp = min (log (p), -1);
  minBF = exp (1 + logp + log (-logp));

  % Calculate the false-positive risk from the minumum Bayes Factor
  L10 = 1 ./ minBF;      % Convert to Maximum Likelihood ratio L10 (P(H1)/P(H0))
  fpr = max (0, 1 ./ (1 + L10));  % Calculate minimum false positive risk 
  fpr(isnan(p)) = NaN; 

end

%--------------------------------------------------------------------------

%% FUNCTION TO PRINT OUTPUT

function print_output (stats, nboot, p, method)

    fprintf (['\nSummary of wild bootstrap null hypothesis significance tests for linear models\n',...
              '*******************************************************************************\n\n']);
    fprintf ('Bootstrap settings: \n');
    fprintf (' Function: pinv (X'' * X) * (X'' * y)\n');
    fprintf (' Resampling method: Wild %sbootstrap-t (H0-imposed)\n', method)
    fprintf (' Number of resamples: %u \n', nboot)
    fprintf (' Standard error calculations:');
    if (isempty (method))
      fprintf (' Heteroscedasticity-Consistent (HC0)\n');
    else
      fprintf (' Cluster Robust (CR0)\n');
    end
    fprintf (' Null value (H0) used for hypothesis testing (p-values): 0 \n')
    fprintf ('\nTest Statistics: \n');
    fprintf (' original    std_err     t-stat     p-val    FPR\n');
    for j = 1:p
      fprintf (' %#-+10.4g  %#-+10.4g  %#-+9.3g', ...
               [stats.original(j), stats.std_err(j), stats.tstat(j)])
      if (stats.pval(j) <= 0.001)
        fprintf ('  <.001');
      elseif (stats.pval(j) < 0.9995)
        fprintf ('   .%03u', round (stats.pval(j) * 1e+03));
      elseif (isnan (stats.pval(j)))
        fprintf ('    NaN');
      else
        fprintf ('   1.000');
      end
      if (stats.fpr(j) <= 0.001)
        fprintf ('  <.001\n');
      elseif (stats.fpr(j) < 0.9995)
        fprintf ('   .%03u\n', round (stats.fpr(j) * 1e+03));
      elseif (isnan (stats.fpr(j)))
        fprintf ('    NaN\n');
      else
        fprintf ('   1.000\n');
      end
    end
    fprintf ('\n');

end

%--------------------------------------------------------------------------

%!demo
%!
%! ## Input univariate dataset
%! heights = [183, 192, 182, 183, 177, 185, 188, 188, 182, 185].';
%!
%! ## Test statistics and p-values (H0 = 0)
%! bootwild(heights);
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input bivariate dataset
%! X = [ones(43,1),...
%!     [01,02,03,04,05,06,07,08,09,10,11,...
%!      12,13,14,15,16,17,18,19,20,21,22,...
%!      23,25,26,27,28,29,30,31,32,33,34,...
%!      35,36,37,38,39,40,41,42,43,44]'];
%! y = [188.0,170.0,189.0,163.0,183.0,171.0,185.0,168.0,173.0,183.0,173.0,...
%!     173.0,175.0,178.0,183.0,192.4,178.0,173.0,174.0,183.0,188.0,180.0,...
%!     168.0,170.0,178.0,182.0,180.0,183.0,178.0,182.0,188.0,175.0,179.0,...
%!     183.0,192.0,182.0,183.0,177.0,185.0,188.0,188.0,182.0,185.0]';
%!
%! ## Compute test statistics and p-values
%! bootwild(y,X);
%!
%! ## Please be patient, the calculations will be completed soon...

%!test
%! ## Test if the mean is equal to population value of 150 (one-tailed test)
%!
%! ## Input univariate dataset
%! H0 = 150;
%! heights = [183, 192, 182, 183, 177, 185, 188, 188, 182, 185].';
%!
%! ## Compute test statistics and p-values
%! stats = bootwild(heights-H0);
%! stats = bootwild(heights-H0,[],[]);
%! stats = bootwild(heights-H0,[],[],2000);
%! stats = bootwild(heights-H0,[],[],2000,1);
%! stats = bootwild(heights-H0,[],[],[],[]);
%! [stats,bootstat] = bootwild(heights);

%!test
%! ## Test if the regression coefficients equal 0
%!
%! ## Input bivariate dataset
%! X = [ones(43,1),...
%!     [01,02,03,04,05,06,07,08,09,10,11,...
%!      12,13,14,15,16,17,18,19,20,21,22,...
%!      23,25,26,27,28,29,30,31,32,33,34,...
%!      35,36,37,38,39,40,41,42,43,44]'];
%! y = [188.0,170.0,189.0,163.0,183.0,171.0,185.0,168.0,173.0,183.0,173.0,...
%!     173.0,175.0,178.0,183.0,192.4,178.0,173.0,174.0,183.0,188.0,180.0,...
%!     168.0,170.0,178.0,182.0,180.0,183.0,178.0,182.0,188.0,175.0,179.0,...
%!     183.0,192.0,182.0,183.0,177.0,185.0,188.0,188.0,182.0,185.0]';
%!
%! ## Compute test statistics and p-values
%! stats = bootwild(y,X);
%! stats = bootwild(y,X,[],2000);
%! stats = bootwild(y,X,[],2000,1);
%! stats = bootwild(y,X,[],[],[]);
%! [stats,bootstat] = bootwild(y,X);