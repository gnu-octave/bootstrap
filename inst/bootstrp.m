%  Function File: bootstrp
%
%  Bootstrap sampling
%
%  bootstat = bootstrp(nboot,bootfun,...)
%  bootstat = bootstrp(...,'Weights',weights)
%  bootstat = bootstrp(...,'Options',paropt)
%  [bootstat,bootsam] = bootstrp(...)
%
%  bootstat = bootstrp(nboot,bootfun,...) draws nboot bootstrap data
%  resamples and returns the statistic computed by bootfun in bootstat
%  [1]. bootfun is a function handle (e.g. specified with @), or a
%  string indicating the function name. The third and later input
%  arguments are data (column vectors, or a matrix), that are used
%  to create inputs for bootfun. The resampling method used throughout
%  is balanced resampling [2].
%
%  bootstat = bootstrp(...,'Weights',weights) specifies observation
%  weights. weights must be a vector of non-negative numbers. The
%  length of weights must be equal to first dimension of the
%  non-scalar input argument(s) to bootfun. The weights are used as
%  bootstrap sampling probabilities. Balanced resampling is extended
%  to resampling with weights [3].
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
%  [bootstat,bootsam] = bootstrp(...) also returns bootidx, a matrix of
%  indices from the bootstrap. Each column in bootsam corresponds
%  to one bootstrap sample and contains the row indices of the values
%  drawn from the nonscalar data to create that sample.
%
%  Bibliography:
%  [1] Efron, and Tibshirani (1993) An Introduction to the
%        Bootstrap. New York, NY: Chapman & Hall
%  [2] Davison et al. (1986) Efficient Bootstrap Simulation.
%        Biometrika, 73: 555-66
%  [3] Booth, Hall and Wood (1993) Balanced Importance Resampling
%        for the Bootstrap. The Annals of Statistics. 21(1):286-298
%
%  bootstrp v2.8.6.0 (29/04/2020)
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
  weights = [];
  paropt = struct;
  paropt.UseParallel = false;
  paropt.nproc = nproc;

  % Assign input arguments to function variables
  nboot = argin1;
  bootfun = argin2;
  argin3 = varargin;
  narg = numel(argin3);
  if narg > 1
    while ischar(argin3{end-1})
      if strcmpi(argin3{end-1},'weights')
        weights = argin3{end};
      elseif any(strcmpi({'options','option'},argin3{end-1}))
        paropt = argin3{end};
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
    [ci, bootstat, S, calcurve, bootsam] = ibootci (nboot, {bootfun, data{:}}, 'Weights', weights,'Options',paropt);
  else
    [ci, bootstat, S, calcurve] = ibootci (nboot, {bootfun, data{:}}, 'Weights', weights,'Options',paropt);
  end
  bootstat = bootstat.';

end

