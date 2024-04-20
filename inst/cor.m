% Vectorized function for computing Pearson's correlation coefficient (RHO)
% between the respective columns in two equal-sized matrices.
%
% -- Function File: RHO = cor (X, Y)
%
%     'RHO = cor (X, Y)' computes Pearson's correlation coefficient (RHO)
%     between the column vectors X and Y. If X and Y are matrices, then
%     RHO will be a row vector corresponding to column-wise correlation
%     coefficients. Hence this function is vectorised for rapid computation
%     of the correlation coefficient in bootstrap resamples. Note that
%     unlike the native @corr function, the correlation coefficients
%     returned here are representative of the finite data sample and are
%     not unbiased estimates of the population parameter.
% 
%     cor (X, Y) = ...
%
%     mean ( (X - mean (X)) .* (Y - mean (Y)) ) ./ (std (X, 1) .* std (Y, 1))
%
%  cor (version 2023.05.02)
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
%  along with this program.  If not, see http://www.gnu.org/licenses/

function RHO = cor (X, Y)

  % Evaluate input arguments
  if ((nargin < 2) || (nargin > 3))
    error ('cor: Invalid number of input arguments')
  end
  if (nargout > 1)
    error ('cor: Invalid number of output arguments')
  end
  szx = size (X);
  szy = size (Y);
  if (~ (szx(1) == szy(1)))
    error ('cor: X and Y must have the same number of rows')
  end
  if (numel (szx) > 2)
    error ('cor: X cannot have more than 2 dimensions')
  end
  if (numel (szy) > 2)
    error ('cor: Y cannot have more than 2 dimensions')
  end

  % Calculate correlation coefficient 
  % Note that the factor of 1/n in the terms for the mean and standard
  % deviations in the numerator and denominator respectively cancel out
  % and so are omitted for efficiency
  XERR = bsxfun (@minus, X, mean (X));
  YERR = bsxfun (@minus, Y, mean (Y));
  if (all (szx == szy))
    RHO = sum ( XERR .* YERR ) ./ sqrt (sum (XERR.^2) .* sum (YERR.^2));
  else
    RHO =  sum (bsxfun (@times, XERR, YERR)) ./ ...
          sqrt (bsxfun (@times, sum (XERR.^2), sum (YERR.^2)));
  end
