%  Function File: bootci
%
%  Bootstrap confidence intervals
%
%  CI = bootci (NBOOT, BOOTFUN, D)
%  CI = bootci (NBOOT, BOOTFUN, D1,...,DN)
%  CI = bootci (NBOOT, {BOOTFUN, D}, Name, Value)
%  CI = bootci (NBOOT, {BOOTFUN, D1,...,DN}, Name, Value)
%  CI = bootci (NBOOT, BOOTFUN, D, Name, Value)
%  CI = bootci (NBOOT, BOOTFUN, D1,...,DN, Name, Value)
%  CI = bootci (...,'type', TYPE)
%  CI = bootci (...,'type', 'stud', 'nbootstd', NBOOTSTD)
%  CI = bootci (...,'type', 'cal', 'nbootcal', NBOOTCAL)
%  CI = bootci (...,'alpha', ALPHA)
%  CI = bootci (...,'seed', SEED)
%  CI = bootci (...,'Options', PAROPT)
%  [CI, BOOTSTAT] = bootci (...)
%
%  CI = bootci (NBOOT, BOOTFUN, D) draws nboot bootstrap resamples from the rows
%  of a data sample D and returns 95% confidence intervals (CI) for the bootstrap 
%  statistics computed by BOOTFUN [1]. BOOTFUN is a function handle (e.g. specified 
%  with @), or a string indicating the function name. The third input argument, 
%  data D (a column vector or a matrix), is used as input for BOOTFUN. The
%  resampling method used throughout is balanced bootknife resampling [2-4].
%
%  CI = bootci (NBOOT, BOOTFUN, D1,...,DN) is as above except that the third and
%  subsequent numeric input arguments are data vectors that are used to create
%  inputs for bootfun.
%
%  CI = bootci (..., 'alpha', ALPHA) where ALPHA sets the lower and upper bounds 
%  of the confidence interval(s). The value of ALPHA must be between 0 and 1.
%  The nominal lower and upper percentiles of the confidence intervals CI are 
%  then 100*(ALPHA/2)% and 100*(1-ALPHA/2)% respectively, and nominal central
%  coverage of the intervals is 100*(1-ALPHA)%. The default value of ALPHA is 0.05.
%
%  CI = bootci (..., 'type', TYPE) computes bootstrap confidence interval CI 
%  using one of the following methods:
%    'per' or 'percentile': Percentile method.
%    'bca': Bias-corrected and accelerated method [5,6] (Default).
%    'stud' or 'student': Studentized (bootstrap-t) confidence interval.
%    'cal': Calibrated percentile method (by double bootstrap [7]).
%  Note that when BOOTFUN is the mean, BCa intervals are automatically expanded
%  using Student's t-distribution in order to improve coverage for small samples
%  [8]. The bootstrap-t method includes an additive correction to stabilize
%  the variance when the sample size is small [9].
%
%  CI = bootci (..., 'type', 'stud', 'nbootstd', NBOOTSTD) computes the
%  Studentized bootstrap bootstrap confidence intervals CI, with the standard
%  error of the bootstrap statistics estimated from NBOOTSTD bootstrap data
%  samples. NBOOTSTD is a positive integer value. The default value of NBOOTSTD
%  is 100.
%
%  CI = bootci (..., 'type', 'cal', 'nbootcal', NBOOTCAL) computes the calibrated
%  percentile bootstrap confidence intervals CI, with the calibrated percentiles
%  of the bootstrap statistics estimated from NBOOTCAL bootstrap data samples.
%  NBOOTCAL is a positive integer value. The default value of NBOOTCAL is 200.
%
%  CI = bootci (..., 'seed', SEED) initialises the Mersenne Twister random number
%  generator using an integer SEED value so that bootci results are reproducible.
%
%  CI = bootci (..., 'Options', PAROPT) specifies options that govern if and how
%  to perform bootstrap iterations using multiple processors (if the Parallel 
%  Computing Toolbox or Octave Parallel package is available). This argument is
%  a structure with the following recognised fields:
%
%   'UseParallel' - If true, use parallel processes to accelerate bootstrap
%                   computations on multicore machines, specifically
%                   non-vectorized function evaluations, double bootstrap
%                   resampling and jackknife function evaluations. Default is
%                   false for serial computation. In MATLAB, the default is
%                   true if a parallel pool has already been started. 
%
%                   Note that the standard error calculations for Studentized
%                   bootstrap confidence intervals are not accelerated by
%                   parallel processing.
%
%   'nproc'       - nproc sets the number of parallel processes
%
%  [CI, BOOTSTAT] = bootci(...) also returns the bootstrap statistics used to
%  calculate the confidence intervals CI.
%
%  [CI, BOOTSTAT, BOOTSAM] = bootci(...) also returns BOOTSAM, a matrix of 
%  indices from the bootstrap. Each column in BOOTSAM corresponds to one 
%  bootstrap sample and contains the row indices of the values drawn from the 
%  nonscalar data argument to create that sample.
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
%  [5] Efron (1987) Better Bootstrap Confidence Intervals. JASA, 
%        82(397): 171-185 
%  [6] Efron, and Tibshirani (1993) An Introduction to the
%        Bootstrap. New York, NY: Chapman & Hall
%  [7] Hall, Lee and Young (2000) Importance of interpolation when
%        constructing double-bootstrap confidence intervals. Journal
%        of the Royal Statistical Society. Series B. 62(3): 479-491
%  [8] Hesterberg, Tim (2014), What Teachers Should Know about the 
%        Bootstrap: Resampling in the Undergraduate Statistics Curriculum, 
%        http://arxiv.org/abs/1411.5279
%  [9] Polansky (2000) Stabilizing bootstrap-t confidence intervals
%        for small samples. Can J Stat. 28(3):501-516
%
%  bootci (version 2022.12.06)
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
  if nargin<2
    error('bootci usage: ''bootci (NBOOT, {BOOTFUN, DATA}, varargin)''; atleast 2 input arguments required');
  end

  % Store local functions in a stucture for parallel processes
  localfunc = struct ('col2args',@col2args,...
                      'empcdf',@empcdf);

  % Check if using MATLAB or Octave
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

  % Apply defaults
  bootfun = @mean;
  alpha = 0.05;
  type = 'bca';
  nbootstd = 100;
  nbootcal = 200;
  paropt = struct;
  paropt.UseParallel = false;
  if ~ISOCTAVE
    ncpus = feature('numcores');
  else
    ncpus = nproc;
  end
  paropt.nproc = ncpus;

  % Assign input arguments to function variables
  nboot = argin1;
  bootfun = argin2;
  argin3 = varargin;
  narg = numel(argin3);
  if narg > 1
    while ischar(argin3{end-1})
      if any(strcmpi({'Options','Option'},argin3{end-1}))
        paropt = argin3{end};
      elseif any(strcmpi('alpha',argin3{end-1}))
        alpha = argin3{end};
      elseif any(strcmpi('type',argin3{end-1}))
        type = argin3{end};
      elseif any(strcmpi('nbootstd',argin3{end-1}))
        nbootstd = argin3{end};
      elseif any(strcmpi('nbootcal',argin3{end-1}))
        nbootcal = argin3{end};
      elseif any(strcmpi('seed',argin3{end-1}))
        seed = argin3{end};
        % Initialise the random number generator with the seed
        boot (1, 1, true, [], seed);
      else
        error('bootci: unrecognised input argument to bootci')
      end
      argin3 = {argin3{1:end-2}};
      narg = numel(argin3);
      if narg < 2
        break
      end
    end
  end
  if iscell(argin2)
    bootfun = argin2{1};
    if (numel(argin2) > 2)
      data = argin2(2:end);
      n = size (data{1},1);
    else
      data = argin2{2};
      n = size (data,1);
    end
  else
    bootfun = argin2;
    if (numel(argin3) > 1)
      data = argin3;
    else
      data = argin3{1};
    end
    n = size (data,1);
  end
  if paropt.UseParallel
    ncpus = paropt.nproc;
  else
    ncpus = 0;
  end

  % Error checking
  if (numel(alpha) > 1)
    error ('bootci: ALPHA must be a scalar value');
  end
  if ~isa (alpha,'numeric')
    error ('bootci: ALPHA must be a numeric');
  end
  if any ((alpha < 0) | (alpha > 1))
    error ('bootci: ALPHA must be a value between 0 and 1');
  end
  if ~isa (nboot, 'numeric')
    error ('bootci: NBOOT must be numeric');
  end
  if (numel (nboot) > 1)
    error ('bootci: NBOOT must be a positive integer');
  end
  if (nboot ~= abs (fix (nboot)))
    error ('bootci: NBOOT must contain positive integers');
  end
  if ~isa (nbootstd, 'numeric')
    error ('bootci: NBOOTSTD must be numeric');
  end
  if (numel (nbootstd) > 1)
    error ('bootci: NBOOTSTD must be a scalar value');
  end
  if (nbootstd ~= abs (fix (nbootstd)))
    error ('bootci: NBOOTSTD must be a positive integer');
  end  
  if ~isa (nbootcal, 'numeric')
    error ('bootci: NBOOTCAL must be numeric');
  end
  if (numel (nbootcal) > 1)
    error ('bootci: NBOOTCAL must be a scalar value');
  end
  if (nbootcal ~= abs (fix (nbootcal)))
    error ('bootci: NBOOTCAL must be a positive integer');
  end    

  % Apply interval type
  switch lower(type)
    case 'bca'
      % Do nothing, BCa intervals are the default in the bootknife function
    case {'per','percentile','stud','student'}
      % Set quantiles directly to calculate percentile intervals
      alpha = [alpha / 2, 1 - alpha / 2];
    case 'cal'
      nboot = cat (2, nboot, nbootcal);
    otherwise
      error ('bootci: interval TYPE not supported')
  end

  % Parse input arguments to the bootknife function to calculate confidence intervals
  switch lower(type)
    case {'stud','student'}
      % Use bootstrap-t method with variance stabilization for small samples
      % Polansky (2000) Can J Stat. 28(3):501-516
      [stats, bootstat, bootsam] = bootknife (data, nboot, bootfun, alpha, [], ncpus);
      % If DATA is a cell array of equal size colunmn vectors, convert the cell
      % array to a matrix and redefine bootfun to parse multiple input arguments
      if iscell(data)
        data = [data{:}];
        bootfun = @(data) localfunc.col2args(bootfun, data);
      end
      % Calculate standard errors of the bootstrap statistics
      bootse = @(BOOTSAM) getfield (bootknife (data(BOOTSAM,:), nbootstd, bootfun, NaN), 'std_error');
      SE = cellfun (bootse, num2cell (bootsam,1));
      a = n^(-3/2) * stats.std_error; % Additive constant to stabilize the variance
      % Calculate Studentized statistics
      ridx = isnan(bootstat); bootstat(ridx) = []; SE(ridx) = [];
      T = (bootstat - stats.original) ./ (SE + a);
      [cdf, T] = localfunc.empcdf (T, 1);
      % Calculate intervals from empirical distribution of the Studentized bootstrap statistics
      tmp = arrayfun ( @(p) stats.original - stats.std_error * interp1 (cdf, T, p, 'linear'), alpha);
      stats.CI_upper = tmp(1); stats.CI_lower = tmp(2);
    otherwise
      [stats, bootstat] = bootknife (data, nboot, bootfun, alpha, [], ncpus);
  end

  % Format output to be consistent with MATLAB's bootci
  ci = [stats.CI_lower; stats.CI_upper];
  bootstat = bootstat.';

end

%--------------------------------------------------------------------------

function retval = col2args (func, x)

  % Usage: retval = col2args (func, x)
  % col2args evaluates func on the columns of x. Each columns of x is passed
  % to func as a separate argument. 

  % Extract columns of the matrix into a cell array
  xcell = num2cell (x, 1);

  % Evaluate column vectors as independent of arguments to bootfun
  retval = func (xcell{:});

end

%--------------------------------------------------------------------------

function [F, x] = empcdf (bootstat, c)

  % Subfunction to calculate empirical cumulative distribution function of bootstat
  %
  % Set c to:
  %  1 to have a complete distribution with F ranging from 0 to 1
  %  0 to avoid duplicate values in x
  %
  % Unlike ecdf, empcdf uses a denominator of N+1

  % Check input argument
  if ~isa(bootstat,'numeric')
    error ('bootknife:empcdf: BOOTSTAT must be numeric');
  end
  if all(size(bootstat)>1)
    error ('bootknife:empcdf: BOOTSTAT must be a vector');
  end
  if size(bootstat,2)>1
    bootstat = bootstat.';
  end

  % Create empirical CDF
  bootstat = sort(bootstat);
  N = sum(~isnan(bootstat));
  [x,F] = unique(bootstat,'last');
  F = F/(N+1);

  % Apply option to complete the CDF
  if c > 0
    x = [x(1);x;x(end)];
    F = [0;F;1];
  end

  % Remove impossible values
  F(isnan(x)) = [];
  x(isnan(x)) = [];
  F(isinf(x)) = [];
  x(isinf(x)) = [];

end

%--------------------------------------------------------------------------

%!demo
%!
%! # Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! # 95% BCa bootstrap confidence intervals for the mean
%! ci = bootci (2000, @mean, data)

%!demo
%!
%! # Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! # 95% calibrated percentile bootstrap confidence intervals for the mean
%! ci = bootci (2000, {@mean, data}, 'type', 'cal','nbootcal',200)
%!
%! # Please be patient, the calculations will be completed soon...

%!demo
%!
%! # Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! # 95% calibrated percentile bootstrap confidence intervals for the median
%! # with smoothing
%! ci = bootci (2000, {@smoothmedian, data}, 'type', 'cal', 'nbootcal', 200)
%!
%! # Please be patient, the calculations will be completed soon...

%!demo
%!
%! # Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! # 90% percentile bootstrap confidence intervals for the variance
%! ci = bootci (2000, {{@var,1}, data}, 'type', 'per', 'alpha', 0.1)

%!demo
%!
%! # Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! # 90% BCa bootstrap confidence intervals for the variance
%! ci = bootci (2000, {{@var,1}, data}, 'type', 'bca', 'alpha', 0.1)

%!demo
%!
%! # Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! # 90% Studentized bootstrap confidence intervals for the variance
%! ci = bootci (2000, {{@var,1}, data}, 'type', 'stud', 'nbootstd', 50, 'alpha', 0.1)
%!
%! # Please be patient, the calculations will be completed soon...

%!demo
%!
%! # Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! # 90% calibrated percentile bootstrap confidence intervals for the variance
%! ci = bootci (2000, {{@var,1}, data}, 'type', 'cal', 'nbootcal', 200, 'alpha', 0.1)
%!
%! # Please be patient, the calculations will be completed soon...

%!demo
%!
%! # Input bivariate dataset
%! x = [2.12,4.35,3.39,2.51,4.04,5.1,3.77,3.35,4.1,3.35, ...
%!      4.15,3.56, 3.39,1.88,2.56,2.96,2.49,3.03,2.66,3]';
%! y  = [2.47,4.61,5.26,3.02,6.36,5.93,3.93,4.09,4.88,3.81, ...
%!       4.74,3.29,5.55,2.82,4.23,3.23,2.56,4.31,4.37,2.4]';
%!
%! # 95% BCa bootstrap confidence intervals for the correlation coefficient
%! ci = bootci (2000, @corr, x, y)
%!
%! # Please be patient, the calculations will be completed soon...

%!test
%! ## Spatial Test Data from Table 14.1 of Efron and Tibshirani (1993)
%! ## An Introduction to the Bootstrap in Monographs on Statistics and Applied 
%! ## Probability 57 (Springer)
%! A = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!      0 33 28 34 4 32 24 47 41 24 26 30 41]';
%! ## Nonparametric 90% percentile confidence intervals (single bootstrap)
%! ## Table 14.2 percentile intervals are 100.8 - 233.9
%! ci = bootci(2000,{{@var,1},A},'alpha',0.1,'type','per','seed',1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   # test boot m-file result
%!   assert (ci(1), 95.32928994082839, 1e-09);
%!   assert (ci(2), 238.4062130177514, 1e-09);
%! end
%! ## Nonparametric 90% BCa confidence intervals (single bootstrap)
%! ## Table 14.2 BCa intervals are 115.8 - 259.6
%! ci = bootci(2000,{{@var,1},A},'alpha',0.1,'type','bca','seed',1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   # test boot m-file result
%!   assert (ci(1), 112.9782684413938, 1e-09);
%!   assert (ci(2), 265.6921865021881, 1e-09);
%! end
%! ## Nonparametric 90% bootstrap-t confidence intervals (double bootstrap)
%! ci = bootci(2000,{{@var,1},A},'alpha',0.1,'type','stud','nbootstd',100,'seed',1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   # test boot m-file result
%!   assert (ci(1), 109.0959028769563, 1e-09);
%!   assert (ci(2), 307.4473656731515, 1e-09);
%! end
%! ## Nonparametric 90% calibrated confidence intervals (double bootstrap)
%! ci = bootci(2000,{{@var,1},A},'alpha',0.1,'type','cal','nbootcal',200,'seed',1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   # test boot m-file result
%!   assert (ci(1), 110.6138073406352, 1e-09);
%!   assert (ci(2), 305.1908284023669, 1e-09);
%! end
%! # Exact intervals based on theory are 118.4 - 305.2 (Table 14.2)
%! # Note that all of the bootknife intervals are slightly wider than the
%! # non-parametric intervals in Table 14.2 because the bootknife (rather than
%! # standard bootstrap) resampling used here reduces small sample bias

%!test
%! ## Data from Table 14.1: Spatial Test Data in DiCiccio and Efron (1996)
%! ## Bootstrap Confidence Intervals. Statistical Science. 11(3):189-228
%! baseline = [2.12,4.35,3.39,2.51,4.04,5.1,3.77,3.35,4.1,3.35, ...
%!             4.15,3.56, 3.39,1.88,2.56,2.96,2.49,3.03,2.66,3]';
%! oneyear  = [2.47,4.61,5.26,3.02,6.36,5.93,3.93,4.09,4.88,3.81, ...
%!             4.74,3.29,5.55,2.82,4.23,3.23,2.56,4.31,4.37,2.4]';
%! ## Nonparametric 90% BCa confidence intervals (single bootstrap)
%! ## Table 2 BCa intervals are 0.55 - 0.85
%! ci = bootci(2000,{@corr,baseline,oneyear},'alpha',0.1,'seed',1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   # test boot m-file result
%!   assert (ci(1), 0.5495318330432346, 1e-09);
%!   assert (ci(2), 0.8460658696851905, 1e-09);
%! end
%! # Exact intervals based on theory are 0.47 - 0.86 (Table 14.2)