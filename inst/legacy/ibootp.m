%  Function File: ibootp
%
%  Bootstrap p-value: One sample test
%
%   p = ibootp(m,bootstat,S,calcurve)
%
%  Performs a two-tailed bootstrap test for the hypothesis that
%  the data used to generate the bootstrap statistics bootstat
%  comes from a distribution with statistic m. m must be a scalar.
%
%  This functon requires the output bootstrap replicate statistics
%  (bootstat) and the output structure (S) from ibootci. Provision
%  of the calibration curve (calcurve) is optional but highly
%  recommended for accurate p-values if the bootstrap method used
%  was the percentile method.
%
%  p-values obtained from ibootp are consistent with the confidence
%  intervals calculated with ibootci.
%
%  When the p-value is too small to calibrate by double bootstrap,
%  this function automatically resorts to calculating the p-value
%  using the Studentized bootstrap (bootstrap-t) with an additive
%  correction factor to stabilize the variance.
%
%  ibootp v1.5.8.0 (05/05/2020)
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


function p = ibootp(m,bootstat,S,calcurve)

  % Check input arguments
  if nargin < 2
    error('ibootp requires atleast 2 input arguments')
  end
  if nargin < 3
    S.type = 'per';
  end

  % Calculate number of bootstrap replicates
  if iscell(bootstat)
    T1 = bootstat{1};
  else
    T1 = bootstat;
  end
  B = numel(T1);

  % Calculate one-sided P value
  switch lower(S.type)
    case {'per','percentile'}
      % Find P-value from the cdf of the bootstrap distribution by linear interpolation
      [cdf,t1] = empcdf(T1,0);
      p = 1-interp1(t1,cdf,m,'linear');
    case {'cper','bca'}
      % Find P-value from the cdf of the bootstrap distribution by linear interpolation
      [cdf,t1] = empcdf(T1,0);
      p = 1-interp1(t1,cdf,m,'linear');
      % Create distribution functions
      stdnormcdf = @(x) 0.5*(1+erf(x/sqrt(2)));
      stdnorminv = @(p) sqrt(2)*erfinv(2*p-1);
      % Bias correction (and acceleration if applicable) to P-value
      z1 = stdnorminv(p);
      p = stdnormcdf(S.z0+((S.z0+z1)/(1-S.a*(S.z0+z1))));
    case {'stud','student'}
      % Use bootstrap-t method
      p = bootstud(m,bootstat,S);
  end

  % Convert to equal-sided two-tailed test
  p = 2*min(p,1-p);

  % Calibration of P-value if applicable
  if nargin > 3 && any(strcmpi(S.type,{'per','percentile','cper','bca'}))
    C = S.nboot(2);
    if C > 0
      idx = sum((calcurve(:,1)<1));
      if 1-p < calcurve(idx,1)
        % Use same calibration of p-value as used for confidence intervals
        calcurve(1,:)=[];calcurve(end,:)=[];
        p = 1 - interp1(calcurve(:,1),calcurve(:,2),1-p,'linear');
      else
        warning(sprintf(['P value is too small to calibrate for this bootstrap distribution. \n'...
          'Switching to p-value calculation by the bootstrap-t method. \n',...
          'Try increasing the number of second bootstrap replicate samples.']));
        % Use bootstrap-t method when p-value is small
        p = bootstud(m,bootstat,S);
        p = 2*min(p,1-p);
      end
    end
  end

  % Check if first bootstrap replicate sample set is large enough
  if (1/p > B/2) || isnan(p)
    warning(sprintf(['P value is too small for this bootstrap distribution. \n'...
            'Try increasing the number of first bootstrap replicate samples.']));
    p = 2/B;
  end

end
