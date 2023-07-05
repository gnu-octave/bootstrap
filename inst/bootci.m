% -- Function File: CI = bootci (NBOOT, BOOTFUN, D)
% -- Function File: CI = bootci (NBOOT, BOOTFUN, D1,...,DN)
% -- Function File: CI = bootci (NBOOT, {BOOTFUN, D}, NAME, VALUE)
% -- Function File: CI = bootci (NBOOT, {BOOTFUN, D1, ..., DN}, NAME, VALUE)
% -- Function File: CI = bootci (...,'type', TYPE)
% -- Function File: CI = bootci (...,'type', 'stud', 'nbootstd', NBOOTSTD)
% -- Function File: CI = bootci (...,'type', 'cal', 'nbootcal', NBOOTCAL)
% -- Function File: CI = bootci (...,'alpha', ALPHA)
% -- Function File: CI = bootci (...,'seed', SEED)
% -- Function File: CI = bootci (...,'Options', PAROPT)
% -- Function File: [CI, BOOTSTAT] = bootci (...)
%
%     'CI = bootci (NBOOT, BOOTFUN, D)' draws nboot bootstrap resamples from
%     the rows of a data sample D and returns 95% confidence intervals (CI) for
%     the bootstrap statistics computed by BOOTFUN [1]. BOOTFUN is a function 
%     handle (e.g. specified with @), or a string indicating the function name. 
%     The third input argument, data D (a column vector or a matrix), is used
%     as input for BOOTFUN. The resampling method used throughout is balanced
%     bootknife resampling [2-4].
%
%     'CI = bootci (NBOOT, BOOTFUN, D1,...,DN)' is as above except that the
%     third and subsequent numeric input arguments are data vectors that are
%     used to create inputs for bootfun.
%
%     'CI = bootci (NBOOT, {BOOTFUN, D}, NAME, VALUE)' is as above but includes
%     setting optional parameters using Name-Value pairs.
%
%     'CI = bootci (NBOOT, {BOOTFUN, D1, ..., DN}, NAME, VALUE)' is as above but
%     includes setting optional parameters using NAME-VALUE pairs.
%
%     bootci can take a number of optional parameters as NAME-VALUE pairs:
%
%     'CI = bootci (..., 'alpha', ALPHA)' where ALPHA sets the lower and upper 
%     bounds of the confidence interval(s). The value of ALPHA must be between
%     0 and 1. The nominal lower and upper percentiles of the confidence
%     intervals CI are then 100*(ALPHA/2)% and 100*(1-ALPHA/2)% respectively,
%     and nominal central coverage of the intervals is 100*(1-ALPHA)%. The
%     default value of ALPHA is 0.05.
%
%     'CI = bootci (..., 'type', TYPE)' computes bootstrap confidence interval 
%     CI using one of the following methods:
%       • 'norm' or 'normal': Using bootstrap bias and standard error [5].
%       • 'per' or 'percentile': Percentile method [1,5].
%       • 'basic': Basic bootstrap method [1,5].
%       • 'bca': Bias-corrected and accelerated method [6,7] (Default).
%       • 'stud' or 'student': Studentized bootstrap (bootstrap-t) [1,5].
%       • 'cal': Calibrated percentile method (by double bootstrap [8]).
%       Note that when BOOTFUN is the mean, BCa intervals are automatically
%       expanded using Student's t-distribution in order to improve coverage
%       for small samples [9]. The bootstrap-t method includes an additive
%       correction to stabilize the variance when the sample size is small [10].
%
%     'CI = bootci (..., 'type', 'stud', 'nbootstd', NBOOTSTD)' computes the
%     Studentized bootstrap confidence intervals CI, with the standard errors
%     of the bootstrap statistics estimated automatically using resampling
%     methods. NBOOTSTD is a positive integer value > 0 defining the number of
%     resamples. Unbiased standard errors are computed using NBOOTSTD bootknife
%     resamples. The default value of NBOOTSTD is 100.
%
%     'CI = bootci (..., 'type', 'cal', 'nbootcal', NBOOTCAL)' computes the
%     calibrated percentile bootstrap confidence intervals CI, with the
%     calibrated percentiles of the bootstrap statistics estimated from NBOOTCAL
%     bootstrap data samples. NBOOTCAL is a positive integer value. The default
%     value of NBOOTCAL is 200.
%
%     'CI = bootci (..., 'seed', SEED)' initialises the Mersenne Twister random
%     number generator using an integer SEED value so that bootci results are
%     reproducible.
%
%     'CI = bootci (..., 'Options', PAROPT)' specifies options that govern if
%     and how to perform bootstrap iterations using multiple processors (if the
%     Parallel Computing Toolbox or Octave Parallel package is available). This
%     argument is a structure with the following recognised fields:
%        • 'UseParallel':  If true, use parallel processes to accelerate
%                          bootstrap computations on multicore machines,
%                          specifically non-vectorized function evaluations,
%                          double bootstrap resampling and jackknife function
%                          evaluations. Default is false for serial computation.
%                          In MATLAB, the default is true if a parallel pool
%                          has already been started. 
%        • 'nproc':        nproc sets the number of parallel processes
%
%     '[CI, BOOTSTAT] = bootci (...)' also returns the bootstrap statistics
%     used to calculate the confidence intervals CI.
%   
%     '[CI, BOOTSTAT, BOOTSAM] = bootci (...)' also returns BOOTSAM, a matrix 
%     of indices from the bootstrap. Each column in BOOTSAM corresponds to one 
%     bootstrap sample and contains the row indices of the values drawn from 
%     the nonscalar data argument to create that sample.
%
%  Bibliography:
%  [1] Efron, and Tibshirani (1993) An Introduction to the
%        Bootstrap. New York, NY: Chapman & Hall
%  [2] Davison et al. (1986) Efficient Bootstrap Simulation.
%        Biometrika, 73: 555-66
%  [3] Booth, Hall and Wood (1993) Balanced Importance Resampling
%        for the Bootstrap. The Annals of Statistics. 21(1):286-298
%  [4] Hesterberg T.C. (2004) Unbiasing the Bootstrap—Bootknife Sampling 
%        vs. Smoothing; Proceedings of the Section on Statistics & the 
%        Environment. Alexandria, VA: American Statistical Association.
%  [5] Davison and Hinkley (1997) Bootstrap Methods and their Application.
%        (Cambridge University Press)
%  [6] Efron (1987) Better Bootstrap Confidence Intervals. JASA, 
%        82(397): 171-185 
%  [7] Efron, and Tibshirani (1993) An Introduction to the
%        Bootstrap. New York, NY: Chapman & Hall
%  [8] Hall, Lee and Young (2000) Importance of interpolation when
%        constructing double-bootstrap confidence intervals. Journal
%        of the Royal Statistical Society. Series B. 62(3): 479-491
%  [9] Hesterberg, Tim (2014), What Teachers Should Know about the 
%        Bootstrap: Resampling in the Undergraduate Statistics Curriculum, 
%        http://arxiv.org/abs/1411.5279
%  [10] Polansky (2000) Stabilizing bootstrap-t confidence intervals
%        for small samples. Can J Stat. 28(3):501-516
%
%  bootci (version 2023.07.04)
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


function [ci, bootstat, bootsam] = bootci (argin1, argin2, varargin)

  % Evaluate the number of function arguments
  if (nargin < 2)
    error ('bootci usage: ''bootci (NBOOT, {BOOTFUN, DATA}, varargin)''; atleast 2 input arguments required');
  end

  % Check if using MATLAB or Octave
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

  % Apply defaults
  alpha = 0.05;
  type = 'bca';
  nbootstd = 100;
  nbootcal = 200;
  paropt = struct;
  paropt.UseParallel = false;
  if (~ ISOCTAVE)
    ncpus = feature ('numcores');
  else
    ncpus = nproc;
  end
  paropt.nproc = ncpus;

  % Assign input arguments to function variables
  nboot = argin1;
  argin3 = varargin;
  narg = numel (argin3);
  if (narg > 1)
    while (ischar (argin3{end-1}))
      if (any (strcmpi ({'Options', 'Option'}, argin3{end-1})))
        paropt = argin3{end};
      elseif (any (strcmpi ('alpha', argin3{end-1})))
        alpha = argin3{end};
      elseif (any (strcmpi ('type', argin3{end-1})))
        type = argin3{end};
      elseif (any (strcmpi ('nbootstd', argin3{end-1})))
        nbootstd = argin3{end};
      elseif (any (strcmpi ('nbootcal', argin3{end-1})))
        nbootcal = argin3{end};
      elseif (any (strcmpi ('seed', argin3{end-1})))
        seed = argin3{end};
        % Initialise the random number generator with the seed
        boot (1, 1, true, seed);
      else
        error ('bootci: Unrecognised input argument to bootci')
      end
      argin3 = {argin3{1:end-2}};
      narg = numel (argin3);
      if (narg < 2)
        break
      end
    end
  end
  if (iscell (argin2))
    bootfun = argin2{1};
    if (numel (argin2) > 2)
      data = argin2(2:end);
      n = size (data{1}, 1);
    else
      data = argin2{2};
      n = size (data, 1);
    end
  else
    bootfun = argin2;
    if (numel (argin3) > 1)
      data = argin3;
      n = size (data{1}, 1);
    else
      data = argin3{1};
      n = size (data, 1);
    end
  end
  if (paropt.UseParallel)
    ncpus = paropt.nproc;
  else
    ncpus = 0;
  end

  % Error checking
  if (numel(alpha) > 1)
    error ('bootci: ALPHA must be a scalar value');
  end
  if (~ isa (alpha,'numeric'))
    error ('bootci: ALPHA must be a numeric');
  end
  if (any ((alpha < 0) | (alpha > 1)))
    error ('bootci: ALPHA must be a value between 0 and 1');
  end
  if (~ isa (nboot, 'numeric'))
    error ('bootci: NBOOT must be numeric');
  end
  if (numel (nboot) > 1)
    error ('bootci: NBOOT must be a positive integer');
  end
  if (nboot ~= abs (fix (nboot)))
    error ('bootci: NBOOT must contain positive integers');
  end
  if (~ isa (nbootstd, 'numeric'))
    error ('bootci: NBOOTSTD must be numeric');
  end
  if (numel (nbootstd) > 1)
    error ('bootci: NBOOTSTD must be a scalar value');
  end
  if (nbootstd < 1)
    error ('bootci: NBOOTSTD must be an integer > 0');
  end  
  if (~ isa (nbootcal, 'numeric'))
    error ('bootci: NBOOTCAL must be numeric');
  end
  if (numel (nbootcal) > 1)
    error ('bootci: NBOOTCAL must be a scalar value');
  end
  if (nbootcal ~= abs (fix (nbootcal)))
    error ('bootci: NBOOTCAL must be a positive integer');
  end
  % If applicable, check we have parallel computing capabilities
  if (ncpus > 1)
    if (ISOCTAVE)
      software = pkg ('list');
      names = cellfun (@(S) S.name, software, 'UniformOutput', false);
      status = cellfun (@(S) S.loaded, software, 'UniformOutput', false);
      index = find (~ cellfun (@isempty, regexpi (names, '^parallel')));
      if (~ isempty (index))
        if (~ logical (status{index}))
          ncpus = 0;
        end
      else
        ncpus = 0;
      end
      if (ncpus == 0)
        % OCTAVE Parallel Computing Package is not installed or loaded
        warning ('bootci:parallel', ...
          'Parallel Computing Package is installed and/or loaded. Falling back to serial processing.');
      end
    else
      info = ver; 
      if (~ ismember ('Parallel Computing Toolbox', {info.Name}))
        ncpus = 0;
      end
      if (ncpus == 0)
        % MATLAB Parallel Computing Toolbox is not installed or loaded
        warning ('bootci:parallel', ...
          'Parallel Computing Toolbox is installed and/or loaded. Falling back to serial processing.');
      end
    end
  end

  % Apply interval type
  switch (lower (type))
    case {'per', 'perc', 'percentile', 'basic', 'norm', 'normal'}
      % Do nothing
    case {'bca', 'stud', 'student'}
      % Set quantiles directly to calculate percentile intervals
      alpha = [alpha / 2, 1 - alpha / 2];
    case 'cal'
      % Set quantiles directly to calibrate confidence intervals
      alpha = [alpha / 2, 1 - alpha / 2];
      nboot = cat (2, nboot, nbootcal);
    otherwise
      error ('bootci: Interval TYPE not supported')
  end

  % Parse input arguments to the bootknife function to calculate confidence intervals
  ci = [];
  switch (lower (type))
 
    case {'norm', 'normal'}
      % Intervals constructed based on bootstrap estimates of bias and standard
      % error assumming that the sampling distribution is Normal
      [stats, bootstat, bootsam] = bootknife (data, nboot, bootfun, NaN, [], ncpus);
      stdnorminv = @(p) sqrt (2) * erfinv (2 * p-1);
      z = stdnorminv (alpha / 2);
      ci = cat (2, stats.original - stats.bias + stats.std_error * z, ...
                   stats.original - stats.bias - stats.std_error * z).';

    case 'basic'
      % The basic bootstrap method.
      % Center bootstrap statistics
      [stats, bootstat, bootsam] = bootknife (data, nboot, bootfun, alpha, [], ncpus);
      ci = cat (2, 2 * stats.original - stats.CI_upper, ...
                   2 * stats.original - stats.CI_lower).';

    case {'stud', 'student'}
      % Use bootstrap-t method with variance stabilization for small samples
      % Polansky (2000) Can J Stat. 28(3):501-516
      [stats, bootstat, bootsam] = bootknife (data, nboot, bootfun, NaN, [], ncpus);
      % Automatically estimate standard errors of the bootstrap statistics
      % Using bootknife resampling
      if (iscell (data))
        % If DATA is a cell array of equal size colunmn vectors, convert the
        % cell array to a matrix and define function to calculate an unbiased 
        % estimate of the standard error using bootknife resampling
        szx = cellfun (@(x) size (x, 2), data);
        data = [data{:}];
        cellfunc = @(bootsam) bootknife (mat2cell (data (bootsam,:), n, szx), nbootstd, bootfun, NaN,  [], 0, [], [], ISOCTAVE);
      else
        cellfunc = @(bootsam) bootknife (data (bootsam,:), nbootstd, bootfun, NaN, [], 0, [], [], ISOCTAVE);
      end
      if (ncpus > 1)
        if (ISOCTAVE)
          % Octave
          % Set unique random seed for each parallel thread
          pararrayfun (ncpus, @boot, 1, 1, false, 1:ncpus);
          bootout = parcellfun (ncpus, cellfunc, num2cell (bootsam, 1), 'UniformOutput', false);
        else
          % MATLAB
          % Set unique random seed for each parallel thread
          parfor i = 1:ncpus; boot (1, 1, false, i); end
          % Perform inner layer of resampling
          bootout = cell (1, nboot(1));
          parfor b = 1:nboot(1); bootout{b} = cellfunc (bootsam(:,b)); end
        end
      else
        bootout = cellfun (cellfunc, num2cell (bootsam, 1), 'UniformOutput', false);
      end
      se = cell2mat (cellfun (@(S) S.std_error, bootout, 'UniformOutput', false));
      % Compute additive constant to stabilize the variance
      a = n^(-3/2) * stats.std_error;
      % Calculate Studentized bootstrap statistics
      T = bsxfun (@minus, bootstat, stats.original) ./ bsxfun (@plus, se, a);
      % Calculate intervals from empirical distribution of the Studentized bootstrap statistics
      m = size (bootstat, 1);
      ci = zeros (m, 2);
      for j = 1:m
        [t, cdf] = bootcdf (T(j,:), true, 1);
        ci(j,1) = stats.original(j) - stats.std_error(j) * interp1 (cdf, t, alpha(2), 'linear', max (t));
        ci(j,2) = stats.original(j) - stats.std_error(j) * interp1 (cdf, t, alpha(1), 'linear', min (t));
      end
      ci = ci.';

    otherwise

      % Other interval types are natively supported in bootknife function
      [stats, bootstat] = bootknife (data, nboot, bootfun, alpha, [], ncpus);

  end

  % Format output to be consistent with MATLAB's bootci
  if (isempty (ci))
    ci = cat (2, stats.CI_lower, stats.CI_upper).';
  end
  bootstat = bootstat.';

end

%--------------------------------------------------------------------------

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 95% BCa bootstrap confidence intervals for the mean
%! ci = bootci (2000, @mean, data)

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 95% calibrated percentile bootstrap confidence intervals for the mean
%! ci = bootci (2000, {@mean, data}, 'type', 'cal','nbootcal',200)
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 95% calibrated percentile bootstrap confidence intervals for the median
%! ## with smoothing
%! ci = bootci (2000, {@smoothmedian, data}, 'type', 'cal', 'nbootcal', 200)
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 90% percentile bootstrap confidence intervals for the variance
%! ci = bootci (2000, {{@var,1}, data}, 'type', 'per', 'alpha', 0.1)

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 90% BCa bootstrap confidence intervals for the variance
%! ci = bootci (2000, {{@var,1}, data}, 'type', 'bca', 'alpha', 0.1)

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! ## 90% Studentized bootstrap confidence intervals for the variance
%! ci = bootci (2000, {{@var,1}, data}, 'type', 'stud', 'nbootstd', 50, 'alpha', 0.1)
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## 90% calibrated percentile bootstrap confidence intervals for the variance
%! ci = bootci (2000, {{@var,1}, data}, 'type', 'cal', 'nbootcal', 200, 'alpha', 0.1)
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Input bivariate dataset
%! x = [2.12,4.35,3.39,2.51,4.04,5.1,3.77,3.35,4.1,3.35, ...
%!      4.15,3.56, 3.39,1.88,2.56,2.96,2.49,3.03,2.66,3].';
%! y  = [2.47,4.61,5.26,3.02,6.36,5.93,3.93,4.09,4.88,3.81, ...
%!       4.74,3.29,5.55,2.82,4.23,3.23,2.56,4.31,4.37,2.4].';
%!
%! ## 95% BCa bootstrap confidence intervals for the correlation coefficient
%! ci = bootci (2000, @cor, x, y)
%!
%! ## Please be patient, the calculations will be completed soon...

%!demo
%!
%! ## Spatial Test Data from Table 14.1 of Efron and Tibshirani (1993)
%! ## An Introduction to the Bootstrap in Monographs on Statistics and Applied 
%! ## Probability 57 (Springer)
%!
%! ## AIM: to construct 90% nonparametric bootstrap confidence intervals for var(A,1)
%! ## var(A,1) = 171.5, and exact intervals based on Normal theory are [118.4, 305.2].
%!
%! ## Calculations using Matlab's 'Statistics and Machine Learning toolbox' (R2020b)
%! ##
%! ## A = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%! ##      0 33 28 34 4 32 24 47 41 24 26 30 41].';
%! ## varfun = @(A) var(A, 1);
%! ## rng('default'); % For reproducibility
%! ## rng('default'); ci1 = bootci (20000,{varfun,A},'alpha',0.1,'type','norm');
%! ## rng('default'); ci2 = bootci (20000,{varfun,A},'alpha',0.1,'type','per');
%! ## rng('default'); ci4 = bootci (20000,{varfun,A},'alpha',0.1,'type','bca');
%! ## rng('default'); ci5 = bootci (20000,{varfun,A},'alpha',0.1,'type','stud');
%! ##
%! ## Summary of results from Matlab's 'Statistics and Machine Learning toolbox' (R2020b)
%! ##
%! ## method             |   0.05 |   0.95 | length | shape |  
%! ## -------------------|--------|--------|--------|-------|
%! ## ci1 - normal       |  108.9 |  247.4 |  138.5 |  1.21 |
%! ## ci2 - percentile   |   97.6 |  235.8 |  138.2 |  0.87 |
%! ## ci4 - BCa          |  114.9 |  260.5 |  145.6 |  1.57 |*
%! ## ci5 - bootstrap-t  |   46.6 |  232.6 |  186.0 |  0.49 |** 
%! ## -------------------|--------|--------|--------|-------|
%! ## parametric - exact |  118.4 |  305.2 |  186.8 |  2.52 |
%! ##
%! ## * Bug in the fx0 subfunction of bootci
%! ## ** Bug in the bootstud subfunction of bootci
%!
%! ## Calculations using the 'boot' and 'bootstrap' packages in R
%! ## 
%! ## library (boot)       # Functions from Davison and Hinkley (1997)
%! ## A <- c(48,36,20,29,42,42,20,42,22,41,45,14,6,0,33,28,34,4,32,24,47,41,24,26,30,41);
%! ## n <- length(A)
%! ##  var.fun <- function (d, i) { 
%! ##                # Function to compute the population variance
%! ##                n <- length (d); 
%! ##                return (var (d[i]) * (n - 1) / n) };
%! ##  boot.fun <- function (d, i) {
%! ##                # Compute the estimate
%! ##                t <- var.fun (d, i);
%! ##                # Compute sampling variance of the estimate using Tukey's jackknife
%! ##                n <- length (d);
%! ##                U <- empinf (data=d[i], statistic=var.fun, type="jack", stype="i");
%! ##                var.t <- sum (U^2 / (n * (n - 1)));
%! ##                return ( c(t, var.t) ) };
%! ## set.seed(1)
%! ## var.boot <- boot (data=A, statistic=boot.fun, R=20000, sim='balanced')
%! ## ci1 <- boot.ci (var.boot, conf=0.90, type="norm")
%! ## ci2 <- boot.ci (var.boot, conf=0.90, type="perc")
%! ## ci3 <- boot.ci (var.boot, conf=0.90, type="basic")
%! ## ci4 <- boot.ci (var.boot, conf=0.90, type="bca")
%! ## ci5 <- boot.ci (var.boot, conf=0.90, type="stud")
%! ##
%! ## library (bootstrap)  # Functions from Efron and Tibshirani (1993)
%! ## set.seed(1); ci4a <- bcanon (A, 20000, var.fun, alpha=c(0.05,0.95))
%! ## set.seed(1); ci5a <- boott (A, var.fun, nboott=20000, nbootsd=500, perc=c(.05,.95))
%! ##
%! ## Summary of results from 'boot' and 'bootstrap' packages in R
%! ##
%! ## method             |   0.05 |   0.95 | length | shape |  
%! ## -------------------|--------|--------|--------|-------|
%! ## ci1  - normal      |  109.4 |  246.8 |  137.4 |  1.21 |
%! ## ci2  - percentile  |   97.8 |  235.6 |  137.8 |  0.87 |
%! ## ci3  - basic       |  107.4 |  245.3 |  137.9 |  1.15 |
%! ## ci4  - BCa         |  116.9 |  264.0 |  147.1 |  1.69 |
%! ## ci4a - BCa         |  116.2 |  264.0 |  147.8 |  1.67 |
%! ## ci5  - bootstrap-t |  111.8 |  291.2 |  179.4 |  2.01 |
%! ## ci5a - bootstrap-t |  112.7 |  292.6 |  179.9 |  2.06 |
%! ## -------------------|--------|--------|--------|-------|
%! ## parametric - exact |  118.4 |  305.2 |  186.8 |  2.52 |
%! 
%! ## Calculations using the 'statistics-bootstrap' package for Octave/Matlab
%! ##
%! ## A = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%! ##      0 33 28 34 4 32 24 47 41 24 26 30 41].';
%! ## ci1 = bootci (20000,{{@var,1},A},'alpha',0.1,'type','norm','seed',1);
%! ## ci2 = bootci (20000,{{@var,1},A},'alpha',0.1,'type','per','seed',1);
%! ## ci3 = bootci (20000,{{@var,1},A},'alpha',0.1,'type','basic','seed',1);
%! ## ci4 = bootci (20000,{{@var,1},A},'alpha',0.1,'type','bca','seed',1);
%! ## ci5 = bootci (20000,{{@var,1},A},'alpha',0.1,'type','stud','nbootstd',100,'seed',1);
%! ## ci6 = bootci (20000,{{@var,1},A},'alpha',0.1,'type','cal','nbootcal',500,'seed',1);
%! ##
%! ## Summary of results from 'statistics-bootstrap' package for Octave/Matlab
%! ##
%! ## method             |   0.05 |   0.95 | length | shape |  
%! ## -------------------|--------|--------|--------|-------|
%! ## ci1 - normal       |  108.2 |  248.5 |  140.3 |  1.22 |
%! ## ci2 - percentile   |   96.6 |  236.7 |  140.1 |  0.87 |
%! ## ci3 - basic        |  106.4 |  246.5 |  140.1 |  1.15 |
%! ## ci4 - BCa          |  115.4 |  266.1 |  150.4 |  1.69 |
%! ## ci5 - bootstrap-t  |  110.8 |  293.5 |  182.8 |  2.01 |
%! ## ci6 - calibrated   |  113.4 |  297.0 |  183.6 |  2.16 |
%! ## -------------------|--------|--------|--------|-------|
%! ## parametric - exact |  118.4 |  305.2 |  186.8 |  2.52 |
%! ##
%! ## Simulation results for constructing 90% confidence intervals for the
%! ## variance of a population N(0,1) from 1000 random samples of size 26
%! ## (analagous to the situation above). Simulation performed using the
%! ## bootsim script with nboot of 2000.
%! ##
%! ## method               | coverage |  lower |  upper | length | shape |
%! ## ---------------------|----------|--------|--------|--------|-------|
%! ## normal               |    84.0% |   3.3% |  12.7% |   0.80 |  1.21 |
%! ## percentile           |    81.1% |   1.2% |  17.7% |   0.77 |  0.92 |
%! ## basic                |    79.5% |   2.6% |  17.9% |   0.78 |  1.09 |
%! ## BCa                  |    85.0% |   4.3% |  10.7% |   0.84 |  1.63 |
%! ## bootstrap-t          |    90.4% |   3.5% |   6.1% |   1.02 |  2.14 |
%! ## calibrated           |    90.7% |   5.1% |   4.2% |   1.13 |  2.73 |
%! ## ---------------------|----------|--------|--------|--------|-------|
%! ## parametric - exact   |    90.8% |   3.7% |   5.5% |   0.99 |  2.52 |
%!
%! ## The equivalent methods for constructing bootstrap intervals in the 'boot'
%! ## and 'bootstrap' packages (in R) and the statistics-bootstrap package (in
%! ## Octave/Matlab) produce intervals with very similar end points, length and
%! ## shape. However, all intervals calculated using the 'statistics-bootstrap'
%! ## package are slightly longer than the intervals calculated in R because
%! ## the 'statistics-bootstrap' package uses bootknife resampling. The scale of
%! ## the sampling distribution for small samples is approximated better by
%! ## bootknife (rather than bootstrap) resampling. 
%!

%!test
%! ## Test for errors when using some different functionalities of bootci
%! warning ('off', 'bootknife:parallel')
%! try
%!   y = randn (20, 1); 
%!   bootci (2000, 'mean', y);
%!   bootci (2000, @mean, y);
%!   bootci (2000, @mean, y, 'alpha', 0.1);
%!   bootci (2000, {'mean', y}, 'alpha', 0.1);
%!   bootci (2000, {@mean, y}, 'alpha', 0.1);
%!   bootci (2000, {@mean, y}, 'alpha', 0.1, 'seed', 1);
%!   bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'norm');
%!   bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'per');
%!   bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'basic');
%!   bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'bca');
%!   bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'stud');
%!   bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'stud', 'nbootstd', 100);
%!   bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'cal');
%!   bootci (2000, {@mean, y}, 'alpha', 0.1, 'type', 'cal', 'nbootcal', 200);
%!   Y = randn (20); 
%!   bootci (2000, 'mean', Y);
%!   bootci (2000, @mean, Y);
%!   bootci (2000, @mean, Y, 'alpha', 0.1);
%!   bootci (2000, {'mean', Y}, 'alpha', 0.1);
%!   bootci (2000, {@mean, Y}, 'alpha', 0.1);
%!   bootci (2000, {@mean, Y}, 'alpha', 0.1, 'seed', 1);
%!   bootci (2000, {@mean, Y}, 'alpha', 0.1, 'type', 'norm');
%!   bootci (2000, {@mean, Y}, 'alpha', 0.1, 'type', 'per');
%!   bootci (2000, {@mean, Y}, 'alpha', 0.1, 'type', 'basic');
%!   bootci (2000, {@mean, Y}, 'alpha', 0.1, 'type', 'bca');
%!   bootci (2000, {@mean, Y}, 'alpha', 0.1, 'type', 'stud');
%!   bootci (2000, {@mean, Y}, 'alpha', 0.1, 'type', 'cal');
%!   y = randn (20,1); x = randn (20,1); X = [ones(20,1),x];
%!   bootci (2000, @cor, x, y);
%!   bootci (2000, @(y,X) pinv(X)*y, y, X);
%!   bootci (2000, @(y,X) pinv(X)*y, y, X, 'alpha', 0.1);
%!   bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1);
%!   bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1, 'type', 'norm');
%!   bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1, 'type', 'per');
%!   bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1, 'type', 'basic');
%!   bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1, 'type', 'bca');
%!   bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1, 'type', 'stud');
%!   bootci (2000, {@(y,X) pinv(X)*y, y, X}, 'alpha', 0.1, 'type', 'cal');
%! catch
%!   warning ('on', 'bootknife:parallel')
%!   rethrow (lasterror)
%! end
%! warning ('on', 'bootknife:parallel')

%!test
%! ## Spatial Test Data from Table 14.1 of Efron and Tibshirani (1993)
%! ## An Introduction to the Bootstrap in Monographs on Statistics and Applied 
%! ## Probability 57 (Springer)
%! A = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!      0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! ## Nonparametric 90% percentile confidence intervals (single bootstrap)
%! ## Table 14.2 percentile intervals are 100.8 - 233.9
%! ci = bootci(2000,{{@var,1},A},'alpha',0.1,'type','per','seed',1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (ci(1), 95.24158837039771, 1e-07);
%!   assert (ci(2), 237.7156378257705, 1e-07);
%! end
%!
%! ## Nonparametric 90% BCa confidence intervals (single bootstrap)
%! ## Table 14.2 BCa intervals are 115.8 - 259.6
%! ci = bootci(2000,{{@var,1},A},'alpha',0.1,'type','bca','seed',1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (ci(1), 113.2388308884533, 1e-07);
%!   assert (ci(2), 264.9901439787903, 1e-07);
%! end
%!
%! ## Nonparametric 90% bootstrap-t confidence intervals (double bootstrap)
%! ci = bootci(2000,{{@var,1},A},'alpha',0.1,'type','stud','nbootstd',100,'seed',1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (ci(1), 109.3907006307414, 1e-07);
%!   assert (ci(2), 306.8096054314711, 1e-07);
%! end
%!
%! ## Nonparametric 90% calibrated percentile confidence intervals (double bootstrap)
%! ci = bootci(2000,{{@var,1},A},'alpha',0.1,'type','cal','nbootcal',200,'seed',1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (ci(1), 110.6138073406352, 1e-07);
%!   assert (ci(2), 305.1908284023669, 1e-07);
%! end
%!
%! ## Exact intervals based on normal theory are 118.4 - 305.2 (Table 14.2)
%! ## Note that all of the bootknife intervals are slightly wider than the
%! ## non-parametric intervals in Table 14.2 because the bootknife (rather than
%! ## standard bootstrap) resampling used here reduces small sample bias

%!test
%! ## Law school data from Table 3.1 of Efron and Tibshirani (1993)
%! ## An Introduction to the Bootstrap in Monographs on Statistics and Applied 
%! ## Probability 57 (Springer)
%! LSAT = [576 635 558 578 666 580 555 661 651 605 653 575 545 572 594].';
%! GPA = [3.39 3.3 2.81 3.03 3.44 3.07 3 3.43 3.36 3.13 3.12 2.74 2.76 2.88 2.96].'; 
%!
%! ## Nonparametric 90% percentile confidence intervals (single bootstrap)
%! ## Percentile intervals on page 266 are 0.524 - 0.928
%! ci = bootci(2000,{@cor,LSAT,GPA},'alpha',0.1,'type','per','seed',1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (ci(1), 0.5056363801008388, 1e-07);
%!   assert (ci(2), 0.9586254199016858, 1e-07);
%! end
%!
%! ## Nonparametric 90% BCa confidence intervals (single bootstrap)
%! ## BCa intervals on page 266 are 0.410 - 0.923
%! ci = bootci(2000,{@cor,LSAT,GPA},'alpha',0.1,'type','bca','seed',1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (ci(1), 0.4119228032301614, 1e-07);
%!   assert (ci(2), 0.9300646701004258, 1e-07);
%! end
%!
%! ## Nonparametric 90% calibrated percentile confidence intervals (double bootstrap)
%! ci = bootci(2000,{@cor,LSAT,GPA},'alpha',0.1,'type','cal','nbootcal',500,'seed',1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   ## test boot m-file result
%!   assert (ci(1), 0.2307337185192847, 1e-07);
%!   assert (ci(2), 0.9444347128107354, 1e-07);
%! end
%! ## Exact intervals based on normal theory are 0.51 - 0.91