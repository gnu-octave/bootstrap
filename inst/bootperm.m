%  Function File: bootperm
%
%  Bootstrap permutation test
%
%  p = bootperm(data,group)
%  p = bootperm(data,group,nboot)
%  p = bootperm(data,group,nboot,bootfun)
%  p = bootperm(data,group,nboot,bootfun,paropt)
%
%  This function provides a bootstrap version of permutation tests for
%  univariate (vector) or multivatiate (matrix) data [1]. The null 
%  hypothesis (H0) is that for a given statistic defined by bootfun, 
%  the data samples (defined by groups) come from the same population. 
%
%  (Note that matrix input for data will mean that calculation of bootfun 
%  will not be vectorized. When this is the case, the function will 
%  return a harmless warning saying so.)
%
%  p = bootperm(data,group) is a k-sample bootstrap permutation test where 
%  group is a vector the same size as y with unique numbers used as group 
%  labels. If the number of groups (k) is 2, this is a 2-sample bootstrap
%  permutation test. If k > 2, this is a bootstrap permutation test for
%  k groups.
%   
%  p = bootperm(data,group,nboot) sets the number of bootstrap resamples. 
%  The default is 5000.
%
%  p = bootperm(data,group,nboot,bootfun) also sets the statistic calculated
%  from the bootstrap samples. This can be a function handle or string
%  corresponding to the function name. The default is @mean or 'mean'.
%
%  p = bootperm(y,...,nboot,bootfun,paropt) specifies options that govern 
%  if and how to perform bootstrap iterations using multiple processors 
%  (if the Parallel Computing Toolbox or Octave Forge package is available).
%  This argument is a structure with the following recognised fields:
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
%   Bibliography:
%   [1] Efron and Tibshirani. Chapter 16 Hypothesis testing with the
%        bootstrap in An introduction to the bootstrap (CRC Press, 1994)
%
%  bootperm v2.0.0.0 (23/12/2021)
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


function p = bootperm(data,group,nboot,bootfun,paropt)

  % Check and process bootperm input arguments
  nvar = size(data,2);
  if (nargin < 2)
    error('bootperm requires atleast two input arguments');
  end
  if size(group,2)>1
    group = group.'; 
  end
  if (numel(group)>1) && (size(data,1) ~= numel(group))
    error('data and group must be vectors the same size')
  end
  if (nargin < 3) || isempty(nboot)
    nboot = 5000;
  end
  if any(size(nboot)>1)
    error('nboot must be scalar. bootperm is not compatible with bootstrap iteration')
  end
  if (nargin < 4) || isempty(bootfun)
    bootfun = 'mean';
  end
  if (nargin < 5)
    paropt = struct;
    paropt.UseParallel = false;
    % Initialise nproc if it doesn't exist
    if ~exist('nproc','builtin')
      nproc = 1;
    end
    paropt.nproc = nproc;
  end

  % Define function to calculate maximum difference among groups
  % H0: Data is exchangeable across all the groups labelled in group
  func = @(data) maxdiff(data,group,bootfun,nvar);
  stat = func(data);

  % Perform resampling and calculate bootstrap statistics
  bootstat = bootstrp(nboot,func,data,'Options',paropt);

  % Calculate p-value
  p = sum(bootstat>stat)/nboot;
  if (p == 0)
    warning('p-value too small to calculate. Try increasing nboot.')
  end

end

%--------------------------------------------------------------------------

function t = maxdiff(Y,g,bootfun,nvar)

  % Calculate maximum difference between bootfun output of all the groups

  % Get size and of the data vector or matrix
  [m,n] = size(Y);
  if (nvar > 1)
    n = 1;
  end

  % Calculate the number (k) of unique groups
  gk = unique(g);
  k = numel(gk);

  % Calculate bootfun on the data from each group
  % (which is simpler and more intuitive than calculating F or MSE)
  Z = zeros(k,n);
  for i = 1:k
    Z(i,:) = feval(bootfun,Y(g==gk(i),:));
  end

  % Calculate maximum difference statistic (t) among the groups
  Z = sort(Z,1);
  t = Z(k,:)-Z(1,:); % sign always positive

end
