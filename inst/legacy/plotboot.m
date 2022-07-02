%  Function File: plotboot
%
%  Plot bootstrap results
%
%   plotboot(bootstat)
%   plotboot(bootstat,ci)
%   plotboot(bootstat,ci,S)
%
%  plotboot v1.0.0.0 (07/06/2020)
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


function plotboot(bootstat,ci,S)

  % Check input arguments
  if nargin > 3
    error('too many input arguments')
  end

  % Calculate number of bootstrap replicates
  if iscell(bootstat)
    T1 = bootstat{1};
  else
    T1 = bootstat;
  end
  B = numel(T1);

  % Plot both Kernel Density Estimate and Histogram
  kdehist(bootstat,'Normal');

  % Plot confidence intervals (if provided)
  if nargin > 1
    ylim=get(gca,'ylim').';
    hold on;
    plot(ones(2,1)*ci(1),ylim,'r:','linewidth',1.5); % Lower CI limit
    plot(ones(2,1)*ci(2),ylim,'r:','linewidth',1.5); % Upper CI limit
    hold off;
  end

  % Plot sample statistic and standard error limits (if ibootci output structure provided)
  if nargin > 2
    hold on;
    plot(ones(2,1)*S.stat,ylim,'k-','linewidth',2.0);      % Sample statistic
    plot(ones(2,1)*S.stat-S.SE,ylim,'r-','linewidth',1.5); % Lower SE limit
    plot(ones(2,1)*S.stat+S.SE,ylim,'r-','linewidth',1.5); % Upper SE limit
    hold off;
  end

end
