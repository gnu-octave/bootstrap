%  Function File: bootanovan
%
%  Bootstrap N-way ANOVA
%
%  p = bootanovan(DATA,GROUP)
%  p = bootanovan(DATA,GROUP)
%  p = bootanovan(DATA,GROUP,nboot)
%  p = bootanovan(DATA,GROUP,nboot,residuals)
%  p = bootanovan(DATA,GROUP,nboot,residuals,nproc)
%  p = bootanovan(DATA,GROUP,nboot,residuals,nproc,varargin)
%  [p,F] = bootanovan(DATA,GROUP,...)
%
%  This bootstrap function is a wrapper for anovan. The p-values are 
%  calculated using bootstrap distributions of the F-statistics.
%  The default approach is to resample the raw data with replacement 
%  assuming exchangeability between the groups across all the factors,
%  akin to Manly's approach of unrestricted permutations [1] but using
%  balanced, bootknife resampling.
%
%  bootanovan requires anovan from either the Statistics package in 
%  Octave or the Statistics and Machine Learning Toolbox in Matlab. 
%  is valid for the sorts of models that Octave anovan can support.
%  Note also that the anovan calculations in Octave do not work unless
%  there is replication for each combination of factor levels. 
%
%  p = bootanovan(DATA,GROUP,nboot) sets the number of bootstrap 
%  resamples. Increasing nboot reduces the Monte Carlo error of the p-
%  value estimates but the calculations take longer to complete. When 
%  nboot is empty or not provided, the default (and minimum allowable 
%  nboot to compute two-tailed p-values down to 0.001) is 1000 - an
%  error is returned if the nboot provided by the user is lower than 
%  this. To reduce monte carlo error and bias, the algorithm uses  
%  balanced, bootknife resampling.
%
%  p = bootanovan(DATA,GROUP,nboot,residuals) sets bootanovan to 
%  resample the ANOVA model residuals (instead of the raw data). 
%  Resampling residuals is akin to the approach of ter Braak for
%  resampling model residuals in permutation and bootstrap tests [1,2].
%  Default is false (i.e. Manly's approach).
%
%  p = bootanovan(DATA,GROUP,nboot,residuals,nproc) sets the number 
%  of parallel processes to use to accelerate computations on 
%  multicore machines. This feature requires the Parallel package 
%  (in Octave), or the Parallel Computing Toolbox (in Matlab).
%
%  p = bootanovan(DATA,GROUP,nboot,residuals,nproc,varargin) allows users to 
%  enter any number of input arguments that will be passed to anovan.
%  These should be key-value pairs of input arguments, for example the 
%  key 'model' supports the following keys:
%   'linear'      - computes N main effects
%   'interaction' - compute N effects and N*(N-1) two-factor interactions
%   'full'        - compute interactions at all levels
%  Please see the help information for anovan for further information 
%  about the key-value pairs supported as input arguments.
% 
%  [p,F] = bootnhst(DATA,GROUP,...) also returns the F-statistics
%
%  [p,F,FDIST] = bootnhst(DATA,GROUP,...) also returns the F-statistics
%  calculated from the bootstrap resamples.
%
%  Bibliography:
%  [1] Howel, D.C. Permutation Tests for Factorial Designs. 
%      Last modified: 03/07/2009, Accessed: 26/07/2022
%      www.uvm.edu/~statdhtx/StatPages/Permutation%20Anova/PermTestsAnova.html
%  [2] ter Braak, CJF (1992) Permutation Versus Bootstrap Significance 
%      Tests in Multiple Regression and ANOVA. In Bootstrapping and Related 
%      Techniques. (K. J. Jockel, Ed.), Springer-Verlag, Berlin, pp. 79-86.
%
%  bootanovan v1.2.0.0 (25/07/2022)
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


function [p, F, FDIST] = bootanovan (data, group, nboot, residuals, ncpus, varargin)

  % Check if running in Octave (else assume Matlab)
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

  % Check for dependency anovan
  if ~exist('anovan','file')
    if ISOCTAVE
      statspackage = ismember ({info.Name}, 'statistics');
      if (~ any (statspackage)) || (str2num(info(statspackage).Version(1:3)) < 1.5)
        error ('bootanovan requires version >= 1.5 of the statistics package')
      end
    else
      error('bootanovan: requires anovan from the Statistics and Machine Learning Toolbox')
    end
  end

  % Check and process bootanovan input arguments
  if (nargin < 2)
    error('bootanovan requires atleast two input arguments');
  end
  if ischar(group)
    group = cellstr(group);
  end
  if (size(group,1)>1) && (size(data,1) ~= size(group,1))
    error('DATA and GROUP must have the same number of rows')
  end
  if iscell(group)
    if ~iscellstr(group)
      group = cell2mat(group);
    end
  end
  if (nargin < 3) || isempty(nboot)
    nboot = 1000;
  else 
    if nboot < 1000
      error('the minimum allowable value of nboot is 1000')
    end
  end
  if any(size(nboot)>1)
    error('nboot must be scalar. bootnhst is not compatible with bootstrap iteration')
  end
  if (nargin < 4)
    residuals = false;
  else
    if ~islogical(residuals) && (numel(residuals) > 1)
      error('residuals must be a logical scalar value')
    end
    if (residuals && ISOCTAVE)
      residuals = false;
      warning('residuals argument is ignored in Octave (residuals = false); the raw data will be resampled.')
    end
  end  
  if (nargin < 5)
    ncpus = 0;    % Ignore parallel processing features
  elseif ~isempty (ncpus) 
    if ~isa (ncpus, 'numeric')
      error('ncpus must be numeric');
    end
    if any (ncpus ~= abs (fix (ncpus)))
      error ('ncpus must be a positive integer')
    end    
    if (numel (ncpus) > 1)
      error ('ncpus must be a scalar value')
    end
  end
  if (nargin < 6)
    options = {};
  else 
    options = varargin;
  end
  if any(strcmpi(options,'alpha'))
    error('the optional anovan parameter ''alpha'' is not supported')
  end
  if (nargout > 3)
    error('bootanovan only supports up to 3 output arguments')
  end

  % Perform balanced, bootknife resampling and compute bootstrap statistics
  boot (1, 1, false, 1, 0); % set random seed to make bootstrap resampling reproducible 
  cellfunc = @(data) anovan_wrapper (data, group, ISOCTAVE, options);
  if residuals
    % Get model residuals, we will resample these instead
    % This is the approach of ter Braak
    [junk1,junk2,stats] = anovan(data,group,'display','off',options{:});
    [junk, FDIST] = bootknife (stats.resid, nboot, cellfunc, [], [], ncpus, [], ISOCTAVE);
    clear junk1 junk2;
  else
    % Resample raw data assuming exchangeability across all factors and factor levels
    % This is the approach of Manly
    [junk, FDIST] = bootknife (data, nboot, cellfunc, [], [], ncpus, [], ISOCTAVE);
  end
  clear junk;

  % Calculate ANOVA F-statistics
  F = cellfunc (data);

  % Calculate p-values
  p = sum (bsxfun (@ge, FDIST, F), 2) / nboot;
  
  % Truncate p-values at the resolution of the test
  res = 1 / nboot;
  p(p<res) = res; 
 
end

