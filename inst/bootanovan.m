%  Function File: bootanovan
%
%  Bootstrap N-way ANOVA
%
%  p = bootanovan(DATA,GROUP)
%  p = bootanovan(DATA,GROUP)
%  p = bootanovan(DATA,GROUP,nboot)
%  p = bootanovan(DATA,GROUP,nboot,varargin)
%  [p,F] = bootanovan(DATA,GROUP,...)
%
%  This bootstrap function is a wrapper for anovan. The p-values are 
%  calculated using from bootstrap distributions of the F-statistics.
%  bootanovan requires anovan from either the Statistics package in 
%  Octave or the Statistics and Machine Learning Toolbox in Matlab. 
%
%  Note that the resampling strategy used depends on whether this
%  function is run in Octave or Matlab. In Octave, the data is 
%  resampled with replacement assuming exchangeability between 
%  the groups across all the factors, which is valid for the sorts 
%  of models that Octave anovan can support. In Matlab, the model 
%  residuals are resampled with replacement making it valid for a  
%  broader range of models (e.g. nesting, random effects, continuos 
%  predictors).
%
%  Note that the anovan calculations in Octave do not work unless
%  there is replication for each combination of factor levels. 
%
%  p = bootanovan(DATA,GROUP,nboot) sets the number of bootstrap 
%  resamples. Increasing nboot reduces the monte carlo error of the p-
%  value estimates but the calculations take longer to complete. When 
%  nboot is empty or not provided, the default (and minimum allowable 
%  nboot to compute two-tailed p-values down to 0.001) is 1000 - an
%  error is returned if the nboot provided by the user is lower than 
%  this. To reduce monte carlo error, the algorithm uses balanced 
%  bootstrap resampling.
%
%  p = bootanovan(DATA,GROUP,nboot,varargin) allows users to 
%  enter any number of input arguments that will be passed to anovan.
%  These should be key-value pairs of input arguments, for example the 
%  key 'model' supports the following keys:
%   'linear'      - computes N main effects
%   'interaction' - compute N effects and N*(N-1) two-factor interactions
%   'full'        - compute interactions at all levels
%  Please see the help information for anovan for further information 
%  about the key-value pairs supported as input arguments. Note that
%  unless specified, the calculations in the ANOVA will use type II
%  sum-of-squares.  
% 
%  [p,F] = bootnhst(DATA,GROUP,...) also returns the F-statistics
%
%  bootanovan v1.1.0.0 (24/01/2022)
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


function [p, F, FDIST] = bootanovan (data, group, nboot, varargin)

  % Check if running in Octave (else assume Matlab)
  info = ver; 
  isoctave = any (ismember ({info.Name}, 'Octave'));

  % Check for dependency anovan
  if ~exist('anovan','file')
    error('missing dependency: anovan')
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
  if nargin < 4
    options = {};
  else 
    options = varargin;
  end
  if ~isoctave
    % Set default sum-of-squares type to II for consistency with Octave
    if ~any(strcmpi(options,'sstype'))
      % Make default sstype 2 if not Octave
      options = cat(2,options,'sstype',2);
    end 
    % Get model residuals
    [junk1,junk2,stats] = anovan(data,group,'display','off',options{:});
    resid = stats.resid;
  end
  if nargout > 3
    error('bootanovan only supports up to 3 output arguments')
  end

  % Perform balanced bootstrap resampling and compute bootstrap statistics
  FDIST = cell(1,nboot);
  m = size(data,1);
  idx = boot (m, nboot, false);
  for h = 1:nboot
    if isoctave
      % OCTAVE: Resample raw data
      FDIST{h} = feval('anovan_wrapper',data(idx(:,h)),group,isoctave,options);
    else
      % MATLAB: Resample model residuals since it is valid for a broader range of models
      FDIST{h} = feval('anovan_wrapper',resid(idx(:,h)),group,isoctave,options);
    end
  end
  FDIST = cell2mat(FDIST);

  % Calculate ANOVA F-statistics
  F = feval('anovan_wrapper',data,group,isoctave,options);

  % Calculate p-values
  p = sum (FDIST >= F (1:end,ones(1,nboot)),2) / nboot;

  % Truncate p-values at the resolution of the test
  res = 1/nboot;
  p(p<res) = res; 
 
end

function  F = anovan_wrapper (y, g, isoctave, options)
  
  if isoctave
    % Octave anovan
    [junk,F] = anovan(y,g,options{:});
  else
    % Matlab anovan
    [junk,tbl] = anovan(y,g,'display','off',options{:});
    F = cell2mat(tbl(2:end,6)); 
  end

end
