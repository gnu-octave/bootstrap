% Function file (boot.m) for generating bootstrap sample indices
%
% bootsam = boot (n, B, u)
%
% INPUT VARIABLES
% n (double) is the number of rows (of the data vector)
% B (double) is the number of bootstrap resamples
% u (boolean) false (for bootstrap) or true (for bootknife)
%
% OUTPUT VARIABLE
% bootsam (uint16) is an n x B matrix of bootstrap resamples
%
% Uniform random numbers are generated using the Mersenne Twister 19937 generator
%
% Author: Andrew Charles Penn (2022)

function bootsam = boot (n, B, u)

  % Error checking
  if (n <= 0) || (n ~= fix(n)) || isinf(n) || isnan(n) || (max (size (n)) > 1)
    error ('n must be a finite positive integer')
  end
  if (n > 2^15-1)
    error ('n exceeds the maximum sample size, 2^15-1')
  end
  if (B <= 0) || (B ~= fix(B)) || isinf(B) || isnan(B) || (max (size (n)) > 1)
    error ('B must be a finite positive integer')
  end
  if (B > realmax('single'))
    error ('B exceeds the maximum number of resamples')
  end
  if (nargin < 3)
    u = 0;
  else
    if ~islogical (u)
      error ('u must be either a false (for bootstrap) or true (for bootknife)')
    end
  end
  
  % Preallocate bootsam
  bootsam = zeros (n, B, 'int16');
  
  % Initialize variable defining the available row counts remaining
  c = ones (n, 1, 'single') * B; 
  
  
  % Perform balanced sampling
  r = 0;
  for b = 1:B
    R = rand (n, 1, 'single');
    if (u)
      % Choose which row of the data to exclude for this bootknife sample
      r = b - fix ((b - 1) / n) * n;
    end
    for i = 1:n
      d = c;  
      if (u)
        d(r) = 0;
      end
      if ~sum (d)
        d = c;
      end
      j = sum (R(i) >= cumsum (d ./sum (d))) + 1;
      bootsam (i, b) = j;
      c(j) = c(j) - 1; 
    end
  end