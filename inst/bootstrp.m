%  Function File: bootstrp
%
%  Bootstrap sampling
%
%  BOOTSTAT = bootstrp (NBOOT, BOOTFUN, D)
%  BOOTSTAT = bootstrp (NBOOT, BOOTFUN, D1,...,DN)
%  BOOTSTAT = bootstrp (..., 'Options', PAROPT)
%  [BOOTSTAT, BOOTSAM] = bootstrp (...) 
%
%  BOOTSTAT = bootstrp (NBOOT, BOOTFUN, D) draws NBOOT bootstrap resamples from
%  the data D and returns the statistic computed by BOOTFUN in BOOTSTAT [1]. 
%  bootstrp resamples from the rows of a data sample D (column vector or a matrix). 
%  BOOTFUN is a function handle (e.g. specified with @), or a string indicating the
%  function name. The third input argument is data (column vector or a matrix),
%  that is used to create inputs for BOOTFUN. The resampling method used 
%  throughout is balanced bootknife resampling [2-4].
%
%  BOOTSTAT = bootstrp (NBOOT, BOOTFUN, D1,...,DN) is as above except that 
%  the third and subsequent numeric input arguments are data vectors that
%  are used to create inputs for BOOTFUN.
%
%  BOOTSTAT = bootstrp (..., 'Options', PAROPT) specifies options that govern if 
%  and how to perform bootstrap iterations using multiple processors (if the 
%  Parallel Computing Toolbox or Octave Parallel package is available). This 
%  argument is a structure with the following recognised fields:
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
%  [BOOTSTAT, BOOTSAM] = bootstrp (...) also returns BOOTSAM, a matrix of indices
%  from the bootstrap. Each column in BOOTSAM corresponds to one bootstrap sample
%  and contains the row indices of the values drawn from the nonscalar data argument 
%  to create that sample.
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
%
%  bootstrp (version 2022.10.08)
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


function [bootstat,bootsam] = bootstrp(argin1,argin2,varargin)

  % Evaluate the number of function arguments
  if nargin<2
    error('bootstrp usage: ''bootstrp (nboot, {bootfun, data}, varargin)''; atleast 2 input arguments required');
  end

  % Check if using MATLAB or Octave
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

  % Apply defaults
  bootfun = @mean;
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
      else
        error('bootstrp: unrecognised input argument to bootstrp')
      end
      argin3 = {argin3{1:end-2}};
      narg = numel(argin3);
      if narg < 3
        break
      end
    end
  end
  bootfun = argin2;
  if (numel(argin3) > 1)
    data = argin3;
  else
    data = argin3{1};
  end
  if ~paropt.UseParallel
    ncpus = 0;
  end


  % Error checking
  if ~all(size(nboot) == [1,1])
    error('bootstrp: nboot must be a scalar value')
  end

  % Parse input arguments to the function bootknife
  [jnk, bootstat, bootsam] = bootknife(data, nboot, bootfun,[],[],nproc);

  % Format output to be consistent with MATLAB's bootstrp
  bootstat = bootstat.';
  bootsam = double(bootsam);

end

