% -- Function File: [x, F] = bootcdf (y)
% -- Function File: [x, F] = bootcdf (y, trim)
% -- Function File: [x, F] = bootcdf (y, trim, m)
% -- Function File: [x, F, P] = bootcdf (...)
%
%     '[x, F] = bootcdf (y)' computes the empirical cumulative distribution
%     function (ECDF) of the vector y of length N. This funtction accounts for
%     the presence of ties and so is suitable for computing the ECDF of
%     bootstrap statistics.
%
%     '[x, F] = bootcdf (y, trim)' removes redundant rows of the ECDF when trim
%     is true. When trim is false, x and F are are the same length as y. The
%     default is true.
%
%     '[x, F] = bootcdf (y, trim, m)' specifies the denominator in the
%     calculation of F as (N + m). Accepted values of m are 0 or 1, with the
%     default being 0. When m is 1, quantiles formed from x and F are akin to
%     qtype 6 in the R quantile function.
%
%     '[x, F, P] = bootcdf (...)' also returns the distribution of P values.
%
%     The syntax in this function code is known to be compatible with
%     recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%     Matlab (v6.5.0 and v7.4.0 on Windows XP).
%
%  bootcdf (version 2023.07.05)
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


function [x, F, P] = bootcdf (y, trim, m)

  % Computes empirical cumulative distribution function and p-value distribution
  % in the presence of ties
  % brainder.org/2012/11/28/competition-ranking-and-empirical-distributions/

  % Check input arguments
  if (nargin > 3)
    error ('bootcdf: too many input arguments provided.');
  end
  if (~ isa (y, 'numeric'))
    error ('bootcdf: y must be numeric.');
  end
  if (all (size (y) > 1))
    error ('bootcdf: y must be a vector.');
  end
  if (size (y, 2) > 1)
    y = y.';
  end
  if (nargin < 2)
    trim = true;
  end
  if ( (~ islogical (trim)) && (~ ismember (trim, [0, 1])) )
    error ('bootcdf: m must be scalar.');
  end
  if (nargin < 3)
    % Denominator in calculation of F is (N + m)
    % When m is 1, quantiles formed from x and F are akin to qtype 6
    % www.rdocumentation.org/packages/stats/versions/3.6.2/topics/quantile
    % Hyndman and Fan (1996) Am Stat. 50(4):361-365
    m = 0;
  end
  if (~ isscalar (m))
    error ('bootcdf: m must be scalar.');
  end
  if (~ ismember (m, [0, 1]))
    error ('bootcdf: m must be either 0 or 1');
  end

  % Check output arguments
  if (nargout > 3)
    error ('bootcdf: too many output arguments requested.');
  end

  % Discard NaN values
  ridx = isnan (y);
  y(ridx) = [];

  % Get size of y
  N = numel (y);

  % Create empirical CDF accounting for ties by competition ranking
  x = sort (y);
  [jnk, IA, IC] = unique (x);
  R = cat (1, IA(2:end) - 1, N);
  F = arrayfun (@(i) R(IC(i)), (1 : N)') / (N + m);

  % Create p-value distribution accounting for ties by competition ranking
  P = 1 - arrayfun (@(i) IA(IC(i)) - 1, (1 : N)') / N;

  % Remove redundancy
  if trim
    M = unique ([x, F, P], 'rows', 'last');
    x = M(:,1); F = M(:,2); P = M(:,3);
  end

end