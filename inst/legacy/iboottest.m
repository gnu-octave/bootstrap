%  Function File: iboottest
%
%  One-sample and paired-sample bootstrap test
%
%   p = iboottest(m,nboot,{bootfun,x})
%   p = iboottest(m,nboot,{bootfun,x},Name,Value)
%   p = iboottest(nboot,{bootfun,x,y},...)
%   [p,ci] = iboottest(...)
%   [p,ci,bootstat] = iboottest(...)
%   [p,ci,bootstat,S] = iboottest(...)
%
%  One-sample or paired-sample bootstrap test for univariate data.
%  The null hypothesis for the paired-sample test is that the
%  bootfun statistic calculated for the difference between x and
%  y is equal to zero. The null hypothesis for the one-sample
%  test is that the bootfun statistic calculated for x is equal
%  to zero. All tests are two-tailed.
%
%  Note that this function resamples the rows of the data x.
%
%  See ibootci documentation for input argument definitions and
%  for Name-Value pairs
%
%  iboottest v1.1.3.0 (01/10/2021)
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


function [p,ci,bootstat,S] = iboottest(arg1,varargin)

  % Check and process iboottest input arguments
  if isa(varargin{1},'cell')
    type = 'paired';
    argin = cat(2,arg1,varargin);
    m = 0;
  else
    type = 'one-sample';
    m = arg1;
    argin = varargin;
  end
  x = argin{2}{2};
  switch type
    case 'one-sample'
      if numel(argin{2})>2
        error('too many data variables provided')
      end
      data = x;
    case 'paired'
      y = argin{2}{3};
      if numel(argin{2})>3
        error('too many data variables provided')
      end
      if any(size(x) ~= size(y))
        error('x and y must have the same size and dimensions')
      end
      data = x-y;
  end
  nboot = argin{1};
  bootfun = argin{2}{1};
  argin(1:2) = [];
  
  % Check number of output arguments requested
  if nargout > 4
    error('Too many output arguments requested')
  end

  % Calculate confidence interval using ibootci
  [ci,bootstat,S,calcurve] = ibootci(nboot,{bootfun,data},argin{:});

  % Calculate p-value using ibootp
  p = ibootp(m,bootstat,S,calcurve);

end
