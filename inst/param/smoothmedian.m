%  Function file: smoothmedian
%
%  Function file: smoothmedian
%
%  [p]  = smoothmedian (x)
%  [p]  = smoothmedian (x, dim)
%  [p]  = smoothmedian (x, dim, Tol)
%  [p, stderr]  = smoothmedian (...)
%
%  If x is a vector, find the univariate smoothed median (p) of x.
%  If x is a matrix, compute the univariate smoothed median value
%  for each column and return them in a row vector.  If the optional
%  argument dim is given, operate along this dimension. Arrays of
%  more than two dimensions are not currently supported. 
%
%  The smoothed median is a slightly smoothed version of the 
%  ordinary median and is an M-estimator that is both robust 
%  efficient:
%
%  | Asymptotic            | Mean |    Median  |    Median  |
%  | properties            |      | (smoothed) | (ordinary) |
%  |-----------------------|------|------------|------------|
%  | Breakdown point       | 0.00 |      0.341 |      0.500 |
%  | Pitman efficacy       | 1.00 |      0.865 |      0.637 |
%
%  Smoothing the median is achieved by minimizing the following
%  objective function:
%
%        S (p) = sum ((x(i) - p).^2 + (x(j) - p).^2) .^ 0.5
%               i < j
% 
%  where i and j refers to the indices of the Cartesian product 
%  of each column of x with itself. No smoothing is carried out
%  when data columns contain one or more NaN or Inf elements.
%  Where this is the case, the ordinary median is returned.
%
%  With the ordinary median as the initial value of p, this function
%  minimizes the above objective function by finding the root of the 
%  first derivative using a fast, but reliable, Newton-Bisection  
%  hybrid algorithm. The tolerance (Tol) is the maximum value of the  
%  first derivative that is acceptable to break from optimization. 
%  The default value of Tol is 1e-03.
%
%  The smoothing works by slightly reducing the breakdown point
%  of the median. Bootstrap confidence intervals using the smoothed
%  median have good coverage for the ordinary median of the
%  population distribution and can be used to obtain second order
%  accurate intervals with Studentized bootstrap and calibrated
%  percentile bootstrap methods [1]. When the population distribution 
%  is thought to be strongly skewed, coverage errors can be reduced 
%  by improving symmetry through appropriate data transformation. 
%  Unlike kernel-based smoothing approaches, bootstrapping smoothmedian 
%  does not require explicit choice of a smoothing parameter or a 
%  probability density function. Since the calculations are accelerated
%  by vectorization, this algorithm has space complexity O(n^2) along 
%  dimension dim, thus it is best-suited for small-to-medium sample
%  sizes, which would benefit most from smoothing.
%  
%  Standard errors of the smoothed median can also be estimated  
%  by requesting the second output argument. Note that the standard 
%  errors are quick and dirty estimates calculated from the first 
%  and second derivatives of the objective function. More reliable  
%  estimates of the standard error can be obtained by bootstrap
%  or bootknife resampling.
%
%  Bibliography:
%  [1] Brown, Hall and Young (2001) The smoothed median and the
%       bootstrap. Biometrika 88(2):519-534
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v6.5.0 and v7.4.0 on Windows XP).
%
%  smoothmedian v1.6.0 (16/12/2021)
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


function [p, stderr] = smoothmedian(x,dim,Tol)

  % Evaluate input arguments
  if (nargin < 1) || (nargin > 3)
    error('Invalid number of input arguments')
  end

  if (nargout > 2)
    error('Invalid number of output arguments')
  end

  if numel(size(x)) > 2
    error('arrays of more than 2 dimensions are not supported')
  end

  % Check data dimensions
  if nargin<2
    if size(x,2)==1
      dim = 1;
    elseif size(x,1)==1
      dim = 2;
    else
      dim = 1;
    end
  end

  if dim<1 || dim>2
    error('dim must be a valid dimension');
  end

  % If applicable, switch dimension
  if dim > 1
    x = x.';
  end

  % Check/set tolerance
  if (nargin < 3) || isempty(Tol)
    Tol = 1e-03;
  else
    if (Tol < 1e-03) && (nargout > 1)
      fprintf('Warning: Setting very low tolerance slows down computations.\n');
    end
  end

  % Check input data type
  if ~isa(x,'double')
    error('the x variable must be double precision')
  end

  % Obtain data dimensions
  s = size(x);
  m = s(1);
  n = s(2);
  l = m*(m-1)/2;

  % Variable transformation to normalize stopping criteria
  % Centre the data on the median and divide by the midrange
  centre = median(x,1);
  for i = find(isnan(centre))
    centre(i) = median(x(~isnan(x(:,i)),i),1);
  end
  midrange = (max(x,[],1)-min(x,[],1))/2;
  x = (x - centre(ones(1,m),:)) ./ (midrange(ones(1,m),:));

  % Find column indices where smoothing is not possible
  no = any(isnan(x)) | any(isinf(x));

  % Obtain m(m-1)/2 pairs from the Cartesian product of each column of
  % x with itself by enforcing the restriction i < j on xi and xj
  q = logical(triu(ones(m,m),1));
  i = uint32((1:m)'*ones(1,m));
  xi = x(i(q),:);
  j = uint32(ones(m,1)*(1:m));
  xj = x(j(q),:);

  % Offset xi and xj by +/- h, where h is small enough to avoid
  % encountering stationary points when calculating function
  % derivatives for the Newton steps.
  h = 4.2819e-06;      % sqrt(0.5*eps^(2/3))
  xi = xi+h;
  xj = xj-h;

  %% Nonlinear root finding by Newton-Bisection hybrid algorithm
  % Set initial bracket bounds
  a = min(x,[],1);
  b = max(x,[],1);
  % Set starting value as the median
  c = zeros(1,n);
  p = c;
  % Initialize and define settings
  idx = 1:n;
  MaxIter = 5000;
  % Calculate commonly used operations and assign them to new variables
  z = (xi-xj).^2;
  y = xi+xj;
  % Start iterations
  for Iter = 1:MaxIter
    % Compute derivatives
    temp = ones(l,1)*c;
    D = (xi-temp).^2+(xj-temp).^2;
    R = sqrt(D);
    T = sum((2*temp-y)./R,1);
    if Iter==1
      T(no) = 0;
    end
    temp = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    % The following calculation of the second derivative(s) is 
    % equivalent to (but much faster to compute than):
    % U = sum ( (xi-xj).^2 .* ((xi-temp).^2 + (xj-temp).^2).^(-3/2) ) 
    U = sum(z.*R./D.^2,1);
    D = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    R = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    % Update bracket bounds
    a(T<-Tol) = c(T<-Tol);
    b(T>+Tol) = c(T>+Tol);
    % Evaluate convergence
    cvg = abs(T)<Tol;
    if any(cvg)
      % Export converged parameters
      p(idx(cvg)) = c(cvg);
      % Avoid excess computations in following iterations
      idx(cvg) = [];
      xi(:,cvg) = [];
      xj(:,cvg) = [];
      z(:,cvg) = [];
      y(:,cvg) = [];
      a(cvg) = [];
      b(cvg) = [];
      c(cvg) = [];
      T(cvg) = [];
      U(cvg) = [];
    elseif all(cvg)
      % Export converged parameters and break from iterations
      p(idx(cvg)) = c(cvg);
      break
    end
    % Compute Newton step (fast quadratic convergence but unreliable)
    W = T./U;
    d = c-W;
    T = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    U = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    W = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    % Compute Bisection step (slow linear convergence but very safe)
    c = 0.5*(a+b);
    % Prefer Newton step if it is within brackets
    nwt = d>a & d<b;
    c(nwt) = d(nwt);
    d = [];  %#ok<NASGU> Reduce memory usage. Faster than using clear.
    nwt =[]; %#ok<NASGU> Reduce memory usage. Faster than using clear.
  end
  if Iter==MaxIter
    fprintf('Warning: Root finding failed to reach the specified tolerance.\n');
    if (nargout > 1)
      fprintf('Warning: Estimates of the standard errors will be more unstable.\n');
    end
  end

  % Backtransform the smoothed median value(s)
  p = p .* midrange + centre;

  % Assign ordinary median where smoothing is not possible
  p(no) = centre(no);

  % Estimate standard error(s) (if requested)
  if nargout > 1
    % Obtain m(m-1)/2 pairs from the Cartesian product of each column of
    % x with itself by enforcing the restriction i < j on xi and xj
    xi = (x(i(q),:) + h) .*... 
             midrange(ones(1,l),:) +...
             centre(ones(1,l),:);   % backtransform data when recreating xi
    xj = (x(j(q),:) - h) .*... 
             midrange(ones(1,l),:) +...
             centre(ones(1,l),:);   % backtransform data when recreating xj
    tmp = p(ones(l,1),:);
    % Calculate standard error(s) for the computed value(s) of the smoothed median
    v0 = (1/(m*(m-1))) * sum((((xi-tmp).^2+(xj-tmp).^2).^(-3/2) .* (xi-xj).^2));
    v  = (1/(m*(m-1))) * sum(((xi+xj-2*tmp)./(((xi-tmp).^2+(xj-tmp).^2).^(0.5))).^2);
    stderr = sqrt((v0.^(-2)) .* v / m);
    % Assign 0 to stderr for x columns with 0 variance 
    stderr(midrange==0) = 0;
  else
    stderr = [];
  end
               
  % If applicable, switch dimension
  if dim > 1
    p = p.';
    stderr = stderr.';
  end
