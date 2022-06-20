%  Function file: smoothmedian
%
%  M  = smoothmedian (x)
%  M  = smoothmedian (x, dim)
%  M  = smoothmedian (x, dim, Tol)
%
%  If x is a vector, find the univariate smoothed median (M) of x.
%  If x is a matrix, compute the univariate smoothed median value
%  for each column and return them in a row vector.  If the optional
%  argument dim is given, operate along this dimension. Arrays of
%  more than two dimensions are not currently supported. 
%
%  The smoothed median is a slightly smoothed version of the 
%  ordinary median and is an M-estimator that is both robust 
%  and efficient:
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
%        S (M) = sum (((x(i) - M).^2 + (x(j) - M).^2).^ 0.5)
%               i < j
% 
%  where i and j refers to the indices of the Cartesian product 
%  of each column of x with itself. No smoothing is carried out
%  when data columns contain one or more NaN or Inf elements.
%
%  With the ordinary median as the initial value of M, this function
%  minimizes the above objective function by finding the root of the 
%  first derivative using a fast, but reliable, Newton-Bisection  
%  hybrid algorithm. The tolerance (Tol) is the maximum value of the  
%  step size that is acceptable to break from optimization. The 
%  default value of Tol = range * 1e-04.
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
%  Bibliography:
%  [1] Brown, Hall and Young (2001) The smoothed median and the
%       bootstrap. Biometrika 88(2):519-534
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v6.5.0 and v7.4.0 on Windows XP).
%
%  smoothmedian v1.8 (19/06/2022)
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


function M = smoothmedian(x,dim,Tol)

  % Evaluate input arguments
  if (nargin < 1) || (nargin > 3)
    error('Invalid number of input arguments')
  end

  if (nargout > 1)
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
  if ~ismember(dim,[1,2])
    error('dim must be a valid dimension');
  end

  % If applicable, switch dimension
  if dim > 1
    x = x.';
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

  % Find column indices where smoothing is not possible
  if any(isnan(x)) | any(isinf(x))
     error("x cannot contain Inf or NaN values")
  end
  
  % Calculate basic statistics for each column of the data
  xmax = max(x,[],1);
  xmin = min(x,[],1);
  range = (xmax - xmin) / 2;
  M = median(x); 
  
  % Check/set tolerance
  if (nargin < 3) || isempty(Tol)
    Tol = range * 1e-04;
  else 
    Tol = Tol * ones(1,n);
  end

  % Obtain m(m-1)/2 pairs from the Cartesian product of each column of
  % x with itself by enforcing the restriction i < j on xi and xj
  q = logical(triu(ones(m,m),1));
  i = uint32((1:m)'*ones(1,m));
  xi = x(i(q),:);
  j = uint32(ones(m,1)*(1:m));
  xj = x(j(q),:);
  idx = 1:n;

  % Nonlinear root finding by Newton-Bisection hybrid algorithm
  % Set starting value as the median
  p = M;
  
  % Set initial bracket bounds
  a = xmin; 
  b = xmax;
  xmin = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
  xmax = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.

  % Calculate commonly used operations and assign them to new variables
  z = (xi-xj).^2;
  y = xi+xj;
  
  % Minimize objective function (vectorized)
  MaxIter = 500;
  for Iter = 1:MaxIter
  
    % Compute first derivative
    temp = ones(l,1)*p;
    D = (xi-temp).^2+(xj-temp).^2;
    R = sqrt(D);
    T = sum((2*temp-y)./R,1);
    temp = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    
    % The following calculation of the second derivative(s) is 
    % equivalent to (but much faster to compute than):
    % U = sum ( (xi-xj).^2 .* ((xi-temp).^2 + (xj-temp).^2).^(-3/2) ) 
    U = sum(z.*R./D.^2,1);
    D = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    R = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    
    % Compute Newton step (fast quadratic convergence but unreliable)
    step = T./U;
    if (Iter==1)
      step(a==b) = 0;
    end
    
    % Evaluate convergence
    cvg = abs(step)<=Tol;
    if any(cvg)
    
      % Export converged parameters
      M(idx(cvg)) = p(cvg);
      
      % Avoid excess computations in following iterations
      idx(cvg) = [];
      xi(:,cvg) = [];
      xj(:,cvg) = [];
      z(:,cvg) = [];
      y(:,cvg) = [];
      a(cvg) = [];
      b(cvg) = [];
      p(cvg) = [];
      step(cvg) = [];
      T(cvg) = [];
      Tol(cvg) = [];
      
    end
    
    % Break from loop when all optimizations have converged
    if all(cvg)
      break
    end
    
    % Update bracket bounds
    a(T<-Tol) = p(T<-Tol);
    b(T>+Tol) = p(T>+Tol);
    
    % Preview new value of the smoothed median
    nwt = p-step;
    
    % Prefer Newton step if it is within brackets
    I = (nwt>a) & (nwt<b);
    p(I) = nwt(I);
    
    % Otherwise, compute Bisection step (slow linear convergence but very safe)
    p(~I) = 0.5 * (a(~I) + b(~I));
    
    % Tidy up
    nwt = [];  %#ok<NASGU> Reduce memory usage. Faster than using clear.
    I = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    T = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    U = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    
  end
  
  if Iter==MaxIter
    fprintf('Warning: Root finding failed to reach the specified tolerance.\n');
    if (nargout > 1)
      fprintf('Warning: Estimates of the standard errors will be more unstable.\n');
    end
  end

  % If applicable, switch dimension
  if dim > 1
    M  = M.';
  end