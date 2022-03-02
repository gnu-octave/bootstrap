%  Function File: bootstrp
%
%  Bootstrap sampling
%
%  bootstat = bootstrp(nboot,bootfun,d)
%  bootstat = bootstrp(nboot,bootfun,d1,...,dN)
%  bootstat = bootstrp(...,'Weights',weights)
%  bootstat = bootstrp(...,'strata',strata)
%  bootstat = bootstrp(...,'cluster',clusters)
%  bootstat = bootstrp(...,'block',blocksize)
%  bootstat = bootstrp(...,'Options',paropt)
%  bootstat = bootstrp(...,'bootsam',bootsam)
%  [bootstat,bootsam] = bootstrp(...)
%
%  bootstat = bootstrp(nboot,bootfun,d,...) draws nboot bootstrap data
%  resamples and returns the statistic computed by bootfun in bootstat
%  [1]. bootfun is a function handle (e.g. specified with @), or a
%  string indicating the function name. The third input argument is data 
%  (column vector or a matrix), that is used to create inputs for bootfun. 
%  The resampling method used throughout is balanced resampling [2].
%
%  bootstat = bootstrp(nboot,bootfun,d1,...,dN) is as above except that 
%  the third and subsequent numeric input arguments are data vectors 
%  that are used to create inputs for bootfun. 
%
%  bootstat = bootstrp(...,'Weights',weights) specifies observation
%  weights. weights must be a vector of non-negative numbers. The
%  length of weights must be equal to first dimension of the
%  non-scalar input argument(s) to bootfun. The weights are used as
%  bootstrap sampling probabilities. Balanced resampling is extended
%  to resampling with weights [3].
%
%  ci = ibootci(nboot,{bootfun,...},...,'strata',strata) specifies a
%  vector containing numeric identifiers of strata. The dimensions of
%  strata must be equal to that of the non-scalar input arguments to
%  bootfun. Bootstrap resampling is stratified so that every stratum is
%  represented in each bootstrap test statistic [4]. If weights are also
%  provided then they are within-stratum weights; the weighting of
%  individual strata depends on their respective sample size.
%
%  ci = ibootci(nboot,{bootfun,...},...,'cluster',clusters) specifies 
%  a column vector (or matrix) of numeric identifiers with the same 
%  number of rows as the data. The identifiers should indicate cluster 
%  membership of the data rows. Whereas strata are fixed, clusters are 
%  resampled. This is achieved by two-stage bootstrap resampling of 
%  residuals with shrinkage correction [4,5,6]. If a matrix is provided 
%  defining additional levels of subsampling in a hierarchical data  
%  model, then level two cluster means are computed and resampled.
%
%  ci = ibootci(nboot,{bootfun,...},...,'block',blocksize) specifies
%  a positive integer defining the block length for block bootstrapping
%  data with serial dependence (e.g. stationary time series). The
%  algorithm uses circular, overlapping blocks. 
%
%  bootstat = bootstrp(...,'Options',paropt) specifies options that 
%  govern if and how to perform bootstrap iterations using multiple 
%  processors (if the Parallel Computing Toolbox or Octave Forge 
%  package is available). This argument is a structure with the 
%  following recognised fields:
%
%   'UseParallel' - If true, compute bootstrap iterations in parallel.
%                   Default is false for serial computation. In MATLAB,
%                   the default is true if a parallel pool has already
%                   been started.
%   'nproc'       - The number of processors to use by Octave. Default
%                   is the number of available processors. If you choose
%                   In Matlab, nproc is ignored and the number of parallel
%                   workers should be predefined beforehand by starting
%                   a parallel pool, else it will use the preferred number
%                   of workers.
%
%  Note that bootstrap resampling is not balanced when operating in 
%  parallel.
%
%  bootstat = bootstrp(...,'bootsam',bootsam) performs bootstrap 
%  computations using the indices from bootsam for the bootstrap 
%  (without the need for further resampling).
%
%  [bootstat,bootsam] = bootstrp(...) also returns bootsam, a  
%  matrix of indices from the bootstrap. Each column in bootsam
%  corresponds to one bootstrap sample and contains the row 
%  indices of the values drawn from the nonscalar data argument 
%  to create that sample.
%
%  Bibliography:
%  [1] Efron, and Tibshirani (1993) An Introduction to the
%        Bootstrap. New York, NY: Chapman & Hall
%  [2] Davison et al. (1986) Efficient Bootstrap Simulation.
%        Biometrika, 73: 555-66
%  [3] Booth, Hall and Wood (1993) Balanced Importance Resampling
%        for the Bootstrap. The Annals of Statistics. 21(1):286-298
%  [4] Davison and Hinkley (1997) Bootstrap Methods and their
%        application. Chapter 3: pg 97-100
%  [5] Gomes et al. (2012) Developing appropriate methods for cost-
%        effectiveness analysis of cluster randomized trials.
%        Medical Decision Making. 32(2): 350-361
%  [6] Ng, Grieve and Carpenter (2013) Two-stage nonparametric
%        bootstrap sampling with shrinkage correction for clustered
%        data. The Stata Journal. 13(1): 141-164
%
%  bootstrp v2.8.7.0 (11/01/2022)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/
%
%  Cite as:
%  Andrew Penn (2019). bootstrp (https://www.github.com/acp29/iboot), GitHub.
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


function [bootstat,bootsam] = bootstrp(argin1,argin2,varargin)

  % Evaluate the number of function arguments
  if nargin<3
    error('Too few input arguments');
  end

  % Initialise nproc if it doesn't exist
  if ~exist('nproc','builtin')
    nproc = 1;
  end

  % Apply defaults
  bootfun = @mean;
  weights = [];
  bootsam = [];
  strata = [];
  clusters = [];
  blocksize = [];
  paropt = struct;
  paropt.UseParallel = false;
  % Initialise nproc if it doesn't exist
  if ~exist('nproc','builtin')
    nproc = 1;
  end
  paropt.nproc = nproc;

  % Assign input arguments to function variables
  nboot = argin1;
  bootfun = argin2;
  argin3 = varargin;
  narg = numel(argin3);
  if narg > 1
    while ischar(argin3{end-1})
      if strcmpi(argin3{end-1},'Weights')
        weights = argin3{end};
      elseif any(strcmpi({'Options','Option'},argin3{end-1}))
        paropt = argin3{end};
      elseif strcmpi(argin3{end-1},'bootsam')
        bootsam = argin3{end};
      elseif any(strcmpi(argin3{end-1},{'strata','stratum','stratified'}))
        strata = argin3{end};
      elseif any(strcmpi(argin3{end-1},{'cluster','clusters'}))
        clusters = argin3{end};
      elseif any(strcmpi(argin3{end-1},{'block','blocks','blocksize'}))
        blocksize = argin3{end};
      else
        error('unrecognised input argument to bootstrp')
      end
      argin3 = {argin3{1:end-2}};
      narg = numel(argin3);
      if narg < 3
        break
      end
    end
  end
  data = argin3;
  
  % Error checking
  if ~all(size(nboot) == [1,1])
    error('nboot must be a scalar value')
  end

  % Parse input arguments to the function ibootci
  try
    pool = gcp('nocreate');
  catch
    pool = [];
  end
  if nargout > 1
    [ci, bootstat, S, calcurve, bootsam] = ibootci (nboot, {bootfun, data{:}},...
                                                    'Weights', weights,...
                                                    'Options',paropt,...
                                                    'bootsam',bootsam,...
                                                    'strata',strata,...
                                                    'cluster',clusters,...
                                                    'block',blocksize);
  else
    [ci, bootstat, S, calcurve] = ibootci (nboot, {bootfun, data{:}},...
                                           'Weights', weights,...
                                           'Options',paropt,...
                                           'bootsam',bootsam,...
                                           'strata',strata,...
                                           'cluster',clusters,...
                                           'block',blocksize);
  end
  bootstat = bootstat.';

end

