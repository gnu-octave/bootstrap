%  Function File: bootci
%
%  Bootstrap confidence intervals
%
%  CI = bootci (NBOOT, BOOTFUN, D)
%  CI = bootci (NBOOT, BOOTFUN, D1,...,DN)
%  CI = bootci (NBOOT,{BOOTFUN, D}, Name,Value)
%  CI = bootci (NBOOT,{BOOTFUN, D1,...,DN},Name,Value)
%  CI = bootci (...,'type', TYPE)
%  CI = bootci (...,'alpha', ALPHA)
%  CI = bootci (...,'Options', PAROPT)
%  [CI, BOOTSTAT] = bootci (...)
%
%  CI = bootci (NBOOT, BOOTFUN, D) draws nboot bootstrap resamples from the rows
%  of a data sample D and returns confidence intervals CI for the bootstrap 
%  statistics computed by BOOTFUN [1]. BOOTFUN is a function handle (e.g. specified 
%  with @), or a string indicating the function name. The third input argument, 
%  data D (a column vector or a matrix), is used as input for BOOTFUN. The
%  resampling method used throughout is balanced bootknife resampling [2-4].
%
%  CI = bootci (NBOOT, BOOTFUN, D1,...,DN) is as above except that the third and
%  subsequent numeric input arguments are data vectors that are used to create
%  inputs for bootfun.
%
%  CI = bootci(..., 'alpha', ALPHA) where ALPHA sets the lower and upper bounds 
%  of the confidence interval(s). The value of ALPHA must be between 0 and 1.
%  The nominal lower and upper percentiles of the confidence intervals CI are 
%  then 100*(ALPHA/2)% and 100*(1-ALPHA/2)% respectively, and nominal central
%  coverage of the intervals is 100*(1-ALPHA)%. 
%
%  CI = bootci(..., 'type', TYPE) computes bootstrap confidence interval CI 
%  using one of the following methods:
%    'per' or 'percentile': Percentile method.
%    'bca': Bias-corrected and accelerated method [5,6] (Default).
%    'cal': Calibrated percentile method (by double bootstrap [7]).
%  Note that when BOOTFUN is the mean, BCa intervals are automatically expanded
%  using Student's t-distribution in order to improve coverage for small samples
%  [8]. 
%
%  CI = bootci(...,'type','cal','nbootcal',NBOOTCAL) computes the calibrated
%  percentile bootstrap confidence intervals CI. The calibrated percentiles of
%  the bootstrap statistics are estimated using bootstrap with NBOOTCAL bootstrap
%  data samples. NBOOTCAL is a positive integer value. The default value of
%  NBOOTCAL is 200.
%
%  CI = bootci(..., 'Options', PAROPT) specifies options that govern if and how
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
%
%  bootci (version 2022.10.08)
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


function [ci,bootstat,bootsam] = bootci(argin1,argin2,varargin)

  % Evaluate the number of function arguments
  if nargin<2
    error('bootci usage: ''bootci (NBOOT, {BOOTFUN, DATA}, varargin)''; atleast 2 input arguments required');
  end

  % Check if using MATLAB or Octave
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

  % Apply defaults
  bootfun = @mean;
  alpha = 0.05;
  type = 'bca';
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
  argin3 = varargin;
  if paropt.UseParallel
    nproc = paropt.nproc;
  else
    nproc = 0;
  end
  if iscell(argin2)
    bootfun = argin2{1};
    if (numel(argin2) > 2)
      data = argin2(2:end);
    else
      data = argin2{2};
    end
    narg = numel(argin3);
    if narg > 1
      while ischar(argin3{end-1})
        if any(strcmpi({'Options','Option'},argin3{end-1}))
          paropt = argin3{end};
        elseif any(strcmpi('alpha',argin3{end-1}))
          alpha = argin3{end};
        elseif any(strcmpi('type',argin3{end-1}))
          type = argin3{end};
        elseif any(strcmpi('nbootcal',argin3{end-1}))
          nbootcal = argin3{end};
        else
          error('bootci: unrecognised input argument to bootci')
        end
        argin3(end-1:end) = [];
        narg = numel(argin3);
        if narg < 2
          break
        end
      end
    end
  else
    bootfun = argin2;
    data = argin3;
  end
  if ~paropt.UseParallel
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
    case {'per','percentile'}
      % Set quantiles directly to calculate percentile intervals
      alpha = [alpha / 2, 1 - alpha / 2];
    case 'cal'
      nboot = cat (2, nboot, nbootcal);
    otherwise
      error ('bootci: interval TYPE not supported')
  end

  % Parse input arguments to the function bootknife
  [stats, bootstat] = bootknife (data, nboot, bootfun, alpha, [], ncpus);

  % Format output to be consistent with MATLAB's bootci
  ci = [stats.CI_lower; stats.CI_upper];
  bootstat = bootstat.';

end

%!test
%! ## Spatial Test Data from Table 14.1 of Efron and Tibshirani (1993)
%! ## An Introduction to the Bootstrap in Monographs on Statistics and Applied 
%! ## Probability 57 (Springer)
%! A = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!      0 33 28 34 4 32 24 47 41 24 26 30 41]';
%! ## Nonparametric 90% percentile confidence intervals (single bootstrap)
%! ## Table 14.2 percentile intervals are 100.8 - 233.9
%! boot (1, 1, true, [], 1); # Set random seed
%! ci = bootci(2000,{{@var,1},A},'alpha',0.1,'type','per');
%! assert (ci(1), 95.19578402366864, 1e-09);
%! assert (ci(2), 238.9609467455621, 1e-09);
%! ## Nonparametric 90% BCa confidence intervals (single bootstrap)
%! ## Table 14.2 BCa intervals are 115.8 - 259.6
%! boot (1, 1, true, [], 1); # Set random seed
%! ci = bootci(2000,{{@var,1},A},'alpha',0.1,'type','bca');
%! assert (ci(1), 115.6455796312253, 1e-09);
%! assert (ci(2), 269.4469269661803, 1e-09);
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
%! boot (1, 1, true, [], 1); # Set random seed
%! ci = bootci(2000,{@corr,baseline,oneyear},'alpha',0.1);
%! assert (ci(1), 0.5477225147834641, 1e-09);
%! assert (ci(2), 0.8457573378934136, 1e-09);
%! # Exact intervals based on theory are 0.47 - 0.86 (Table 14.2)