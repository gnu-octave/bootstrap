% Function file (boot.m) for generating bootstrap sample indices
%
% bootsam = boot (n, nboot, u)
%
% INPUT VARIABLES
% n is the number of rows (of the data vector)
% nboot is the number of bootstrap resamples
% u (boolean) for unbiased: false (for bootstrap) or true (for bootknife)
%
% OUTPUT VARIABLE
% bootsam (uint16) is an n x nboot matrix of bootstrap resamples
%
% Uniform random numbers are generated using the Mersenne Twister 19937 generator
%
% Author: Andrew Charles Penn (2022)

function bootsam = boot (n, nboot, u)

  % Error checking
  if (n <= 0) || (n ~= fix(n)) || isinf(n) || isnan(n) || (max (size (n)) > 1)
    error ('n must be a finite positive integer')
  end
  if (n > 2^15-1)
    error ('n exceeds the maximum sample size, 2^15-1')
  end
  if (nboot <= 0) || (nboot ~= fix(nboot)) || isinf(nboot) || isnan(nboot) || (max (size (n)) > 1)
    error ('nboot must be a finite positive integer')
  end
  if (nboot > realmax('single'))
    error ('nboot exceeds the maximum number of resamples')
  end
  if (nargin < 3)
    u = 0;
  else
    if ~islogical (u)
      error ('u must be either a false (for bootstrap) or true (for bootknife)')
    end
  end
  
  % Preallocate bootsam
  bootsam = zeros (n, nboot, 'int16');
  
  % Initialize variable defining the available row counts remaining
  c = ones (n, 1, 'single') * nboot; 
  
  % Perform balanced sampling
  r = 0;
  for b = 1:nboot
    R = rand (n, 1, 'single');
    if (u)
      % Choose which row of the data to exclude for this bootknife sample
      if (fix ((b - 1) / n) == fix (nboot / n))
        r = 1 + fix (rand (1, 1, 'single') * n);  % random
      else
        r = b - fix ((b - 1) / n) * n;            % systematic
      end
    end
    for i = 1:n
      d = c;  
      if (u)
        d(r) = 0;
      end
      if ~sum (d)
        d = c;
      end
      d = cumsum (d);
      j = sum (R(i) >= d ./ d(end)) + 1;
      bootsam (i, b) = j;
      c(j) = c(j) - 1; 
    end
  end