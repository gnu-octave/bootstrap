% Function file for generating bootstrap sample indices
%
% bootsam = boot (n, B, w)
%
% INPUT VARIABLES
% n (double) is the number of rows (of data vector)
% B (double) is the number of bootstrap resamples
% w (double) should be set to 0 (for bootstrap) or 1 (for bootknife)
%
% OUTPUT VARIABLE
% bootsam (uint16) is an n x B matrix of bootstrap resamples
%
% Author: Andrew Charles Penn (2022)

function bootsam = boot (n, B, w)

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
    w = 0;
  else
    if ~ismember (w, [0, 1]) || (max (size (n)) > 1)
     error ('w must be either a 0 (for bootstrap) or a 1 (for bootknife)')
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
    if (w > 0)
      % Choose which row of the data to exclude for this bootknife sample
      r = b - fix ((b - 1) / n) * n;
    end
    for i = 1:n
      d = c;  
      if (w > 0)
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