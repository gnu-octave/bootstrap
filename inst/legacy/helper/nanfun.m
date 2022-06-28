function Y = nanfun (func, X, dim)

  % Helper function file required for ibootci

  % Math functions, ignoring NaNs.
  if nargin < 3
    dim = 1;
  end
  [m,n] = size(X);
  if dim == 1
    [m,n] = size(X);
    Y = zeros(1,n);
  elseif dim == 2
    [n,m] = size(X);
    Y = zeros(m,1);
  else 
    error('dim input argument should be 1 or 2')
  end
  if all(size(X)>1)
    for i = 1:n
      if dim == 1
        Y(1,i) = feval(func, X(~isnan(X(:,i)),i));
      elseif dim == 2
        Y(i,1) = feval(func, X(i,~isnan(X(i,:))));
      end
    end
  else
    Y = feval(func, X(~isnan(X)));
  end

end
