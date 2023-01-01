%  Function file: cor
%
%  RHO = cor (X, Y, TYPE)
%
%  If X and Y are column vectors, a scalar value is returned represnting the 
%  correlation coefficient RHO. If X and Y are matrices, then RHO will be a
%  row vector corresponding to column-wise correlation coefficients. Hence this
%  function is vectorised for rapid computation of the correlation coefficient
%  in bootstrap resamples. Note that unlike the native @corr function, the
%  correlation coefficients returned here are representative of the finite
%  data sample and are not unbiased estimates of the population parameter.
% 
%  cor (X, Y) = ...
%
%   mean ( (X - mean (X)) .* (Y - mean (Y)) ) ./ (std (X, 1) .* std (Y, 1))
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v6.5.0 and v7.4.0 on Windows XP).
%
%  cor (version 2023.01.01)
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


function RHO = cor (X, Y, TYPE)

  % Evaluate input arguments
  if ((nargin < 2) || (nargin > 3))
    error ('cor: Invalid number of input arguments')
  end
  if (nargout > 1)
    error ('cor: Invalid number of output arguments')
  end
  sz = size (X); n = sz (1);
  if (~ all (sz == size (Y)))
    error ('cor: X and Y must be the same size')
  end
  if (numel (sz) > 2)
    error ('cor: arrays of more than 2 dimensions are not supported')
  end

  % If TYPE is Spearman, convert X and Y to ranks
  if (nargin < 3)
    TYPE = 'Pearson';
  end
  switch (lower (TYPE))
  case {'spearman', 'spearmans'}
    for i = 1:sz(2)
      X(:,i) = tiedrank (X(:,i));
      Y(:,i) = tiedrank (Y(:,i));
    end
  case {'pearson', 'pearsons'}
    % Do nothing
  end

  % Calculate correlation coefficient
  XERR = X - mean (X);
  YERR = Y - mean (Y);
  RHO = sum ( XERR .* YERR ) ./ (sqrt (sum (XERR.^2)) .* sqrt (sum (YERR.^2) ));
