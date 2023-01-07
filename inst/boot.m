% Function file for generating balanced bootstrap sample indices or for 
% generating balanced bootstrap resamples of a data vector
%
% USAGE
% BOOTSAM = boot (N, NBOOT)
% BOOTSAM = boot (X, NBOOT)
% BOOTSAM = boot (..., NBOOT, UNBIASED)
% BOOTSAM = boot (..., NBOOT, UNBIASED, SEED)
% BOOTSAM = boot (..., NBOOT, UNBIASED, SEED, WEIGHTS)
%
% INPUT VARIABLES
% N (double) is the number of rows (of the data vector)
% X (double) is a data vector intended for resampling
% NBOOT (double) is the number of bootstrap resamples
% UNBIASED (boolean): false (for bootstrap) or true (for bootknife)
% SEED (double) is a seed for the pseudo-random number generator. 
% WEIGHTS (double) is a weight vector of length n. 
%
% OUTPUT VARIABLE
% bootsam (double) is an n x nboot matrix of resampled data or indices
%
% NOTES
% UNBIASED is an optional input argument. The default is false. If UNBIASED is
% true, then the sample index for omission in each bootknife resample is
% selected systematically. When the remaining number of bootknife resamples is
% not divisible by the sample size (N), then the sample index omitted is
% selected randomly. 
% SEED is an optional scalar input argument used to initialize the random
% number generator to make resampling reproducible between calls to boot.
% WEIGHTS is an optional input argument. If WEIGHTS is empty or not provided,
% the default is a vector of each element equal to NBOOT (i.e. uniform
% weighting). Each element of WEIGHTS is the number of times that the
% corresponding index is represented in bootsam. Therefore, the sum of WEIGHTS
% should equal N * NBOOT. 
%
% Note that the mex function compiled from this source code is not thread safe.
% Below is an example of a line of code one can run in Octave/Matlab before
% attempting parallel operation of boot.mex in order to ensure that the initial
% random seeds of each thread are unique:
%
% In Octave:
% >> pararrayfun(nproc, @boot, 1, 1, false, 1:nproc)
%
% In Matlab:
% >> ncpus = feature('numcores'); parfor i = 1:ncpus; boot (1, 1, false, i); end;
%
% Author: Andrew Charles Penn (2022)


function bootsam = boot (x, nboot, u, s, w)

  % Input variables
  n = numel(x);
  if (n > 1)
    sz = size(x);
    isvec = true;
    if all(sz > 1)
      error('The first input argument must be either a scalar (N) or vector (X).');
    end
  else
    n = x;
    isvec = false;
    if ( (n <= 0) || (n ~= fix(n)) || isinf(n) || isnan(n) )
      error ('The first input argument must be a finite positive integer.')
    end
  end
  if (nboot <= 0) || (nboot ~= fix(nboot)) || isinf(nboot) || isnan(nboot) || (max (size (nboot)) > 1)
    error ('The second input argument (NBOOT) must be a finite positive integer')
  end
  if (nargin > 2) && ~isempty(u)
    if ~islogical (u)
      error ('The third input argument (UNBIASED) must be a logical scalar value')
    end
  else
    u = false;
  end
  if (nargin > 3) && ~isempty(s)
    if (isinf(s) || isnan(s) || (max (size (s)) > 1))
      error ('The fourth input argument (SEED) must be a finite scalar value')
    end
    rand ('twister', s);
  end

  % Preallocate bootsam
  bootsam = zeros (n, nboot);

  % Initialize weight vector defining the available row counts remaining
  if (nargin > 4) && ~isempty(w)
    % Assign user defined weights (counts)
    % Error checking
    if (numel(w) ~= n)
      error('WEIGHTS must be a vector of length N or be the same length as X.');
    end
    if (sum(w) ~= n * nboot)
      error('The elements of WEIGHTS must sum to N * NBOOT.')
    end
    c = w;
  else
    % Assign weights (counts) for uniform sampling
    c = ones (n, 1) * nboot; 
  end

  % Perform balanced sampling
  r = 0;
  for b = 1:nboot
    R = rand (n, 1);
    if (u)
      % Choose which row of the data to exclude for this bootknife sample
      if (fix ((b - 1) / n) == fix (nboot / n))
        r = 1 + fix (rand (1) * n);     % random
      else
        r = b - fix ((b - 1) / n) * n;  % systematic
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
      if (isvec) 
        bootsam (i, b) = x(j);
      else
        bootsam (i, b) = j;
      end
      c(j) = c(j) - 1; 
    end
  end


%!demo
%!
%! % N as input; balanced bootstrap resampling with replacement
%! boot(3, 20, false)

%!demo
%!
%! % N as input; (unbiased) balanced bootknife resampling with replacement
%! boot(3, 20, true)

%!demo
%! 
%! % N as input; balanced resampling with replacement; setting the random seed
%! boot(3, 20, false, 1) % Set random seed
%! boot(3, 20, false, 1) % Reset random seed, BOOTSAM is the same
%! boot(3, 20, false)    % Without setting random seed, BOOTSAM is different

%!demo
%! % Vector (X) as input; balanced resampling with replacement; setting weights
%! x = [23; 44; 36];
%! boot(x, 10, false, 1)            % equal weighting
%! boot(x, 10, false, 1, [20;0;10]) % unequal weighting, no x(2) in BOOTSAM 

%!test
%! ## Test that random samples vary between calls to boot.
%! I1 = boot (3, 20);
%! I2 = boot (3, 20);
%! assert (all (I1(:) == I2(:)), false);

%!test
%! ## Test that random seed gives identical resamples when UNBIASED is false.
%! I1 = boot (3, 20, false, 1);
%! I2 = boot (3, 20, false, 1);
%! assert (all (I1(:) == I2(:)), true);

%!test
%! ## Test that random seed gives identical resamples when UNBIASED is true.
%! I1 = boot (3, 20, true, 1);
%! I2 = boot (3, 20, true, 1);
%! assert (all (I1(:) == I2(:)), true);

%!test
%! ## Test that default setting for UNBIASED is false.
%! I1 = boot (3, 20, [], 1);
%! I2 = boot (3, 20, false, 1);
%! assert (all (I1(:) == I2(:)), true);

%!test
%! ## Test that resampling is balanced when UNBIASED is false.
%! I = boot (3, 20, false, 1);
%! assert (sum (I(:) == 1), 20, 1e-03);
%! assert (sum (I(:) == 2), 20, 1e-03);
%! assert (sum (I(:) == 3), 20, 1e-03);

%!test
%! ## Test that resampling is balanced when UNBIASED is true.
%! I = boot (3, 20, true, 1);
%! assert (sum (I(:) == 1), 20, 1e-03);
%! assert (sum (I(:) == 2), 20, 1e-03);
%! assert (sum (I(:) == 3), 20, 1e-03);

%! ## Test for unbiased sampling (except for last sample).
%! ## The exception is a requirement to maintain balance.
%! I = boot (3, 20, true, 1);
%! assert (all (diff (sort (I(1:end-1)))), false);

%!test
%! ## Test feature for changing resampling weights when UNBIASED is false
%! I = boot (3, 20, false, 1, [30,30,0]);
%! assert (any (I(:) == 3), false);

%!test
%! ## Test feature for changing resampling weights when UNBIASED is true
%! I = boot (3, 20, true, 1, [30,30,0]);
%! assert (any (I(:) == 3), false);