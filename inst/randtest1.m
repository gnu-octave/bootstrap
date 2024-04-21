% Performs a permutation or randomization test to assess if a sample comes from
% a population with a value given for the mean or some other location parameter 
%
% -- Function File: PVAL = randtest1 (A, M)
% -- Function File: PVAL = randtest1 (A, M, NREPS)
% -- Function File: PVAL = randtest1 (A, M, NREPS)
% -- Function File: PVAL = randtest1 (A, M, NREPS, FUNC)
% -- Function File: PVAL = randtest1 (A, M, NREPS, FUNC, SEED)
% -- Function File: PVAL = randtest1 ([A, GA], ...)
% -- Function File: [PVAL, STAT] = randtest1 (...)
% -- Function File: [PVAL, STAT, FPR] = randtest1 (...)
% -- Function File: [PVAL, STAT, FPR, PERMSTAT] = randtest1 (...)
%
%     'PVAL = randtest1 (A, M)' performs a randomization (or permutation) test
%     to ascertain whether data sample in the column vector A comes from a
%     population with mean equal to the value M. The value returned is a 2-
%     tailed p-value against the null hypothesis computed using the absolute
%     values of the mean. This function generates resamples by independently
%     and randomly flipping the signs of values in (A - M).
%
%     'PVAL = randtest1 (A, M, NREPS)' specifies the number of resamples to
%     take in the randomization test. By default, NREPS is 5000. If the number
%     of possible permutations is smaller than NREPS, the test becomes exact.
%     For example, if the number of sampling units (i.e. rows) in the sample
%     is 12, then the number of possible permutations is 2^12 = 4096, so NREPS
%     will be truncated at 4096 and sampling will systematically evaluate all
%     possible permutations. 
%
%     'PVAL = randtest1 (A, M, NREPS, FUNC)' specifies a custom function
%     calculated on the original samples, and the permuted or randomized
%     resamples. Note that FUNC must compute a location parameter and
%     should either be a:
%        o function handle or anonymous function,
%        o string of function name, or
%        o a cell array where the first cell is one of the above function
%          definitions and the remaining cells are (additional) input arguments 
%          to that function (other than the data arguments).
%        See the built-in demos for example usage using the mean.
%
%     'PVAL = randtest1 (A, M, NREPS, FUNC, SEED)' initialises the Mersenne
%     Twister random number generator using an integer SEED value so that
%     the results of 'randtest1' are reproducible when the test is approximate
%     (i.e. when using randomization if not all permutations can be 
%     evaluated systematically).
%
%     'PVAL = randtest1 ([A, GA], M, ...)' also specifies the sampling
%     units (i.e. clusters) using consecutive positive integers in GA for A.
%     Defining the sampling units has applications for clustered resampling,
%     for example in the cases of nested experimental designs. Note that when
%     sampling units contain different numbers of values, function evaluations
%     after sampling cannot be vectorized. If the parallel computing toolbox
%     (Matlab) or parallel package (Octave) is installed and loaded, then the
%     function evaluations will be automatically accelerated by parallel
%     processing on platforms with multiple processors.
%
%     '[PVAL, STAT] = randtest1 (...)' also returns the test statistic.
%
%     '[PVAL, STAT, FPR] = randtest1 (...)' also returns the minimum false
%     positive risk (FPR) calculated for the p-value, computed using the
%     Sellke-Berger approach.
%
%     '[PVAL, STAT, FPR, PERMSTAT] = randtest1 (...)' also returns the
%     statistics of the permutation distribution.
%
%  randtest1 (version 2024.04.21)
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

function [pval, stat, fpr, STATS] = randtest1 (x, m, nreps, func, seed)


  % Check if we are running Octave or Matlab
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

  % Check the number of function arguments
  if (nargin < 2)
    error ('randtest1: A and m must be provided')
  end
  if (nargin > 5)
    error ('randtest1: Too many input arguments')
  end
  if (nargout > 4)
    error ('randtest1: Too many output arguments')
  end

  % Set defaults
  if ( (nargin < 3) || isempty (nreps) )
    nreps = 5000;
  end
  if ( (nargin > 3) && (~ isempty (func)) )
    if (iscell (func))
      args = func(2:end);
      if (ischar (func{1}))
        % Convert character string of a function name to a function handle
        func = str2func (func{1});
      else
        func = func{1};
      end
    elseif (ischar (func))
      % Convert character string of a function name to a function handle
      func = str2func (func);
      args = {};
    elseif (isa (func, 'function_handle'))
      args = {};
    else
      error ('randtest1: FUNC must be a function name or function handle')
    end
  else
    func = @mean;
    args = {};
  end
  if ( (nargin > 4) && (~ isempty (seed)) )
    % Set random seed
    rand ('seed', seed);
  end
  func2 = @(A, B) func (A, args{:});  % func in randtest2 takes two data
                                      % arguments but we only want to evaluate
                                      % the function on one of them

  % Error checking
  if ((~ isscalar (m)) || isinf (m) || isnan (m) || ...
      (~ or (isa (m, 'single'), isa (m, 'double'))))
    error (cat (2, 'randtest1: The second input argument (M) must be a', ...
                   ' finite scalar value.'))
  end
  
  % Get size of the data
  szx = size (x);
  if (numel (szx) > 2)
    error ('randtest1: Variable A has too many dimensions')
  end
  
  % Subtract the reference value from the data
  if (szx(2) > 1)
    % Case when clusters defined in second
    gx = x(:, 2);
    z  = cat (2, x(:, 1) - m, gx);
    zf = cat (2, m - x(:, 1), gx);
  else
    z = x - m;
    zf = -z;
  end

  % Perform ranomization or permutation test
  % We can achieve this using randtest2 with paired == true since a one-sample
  % randomization/permutation test is a special case of a paired two-sample
  % test designed to flip the sign of the sample values
  [pval, stat, fpr, STATS] = randtest2 (z, zf, true, nreps, func2);

end

%!demo
%!
%! % Mouse data from Table 2 (page 11) of Efron and Tibshirani (1993)
%! treatment = [94 197 16 38 99 141 23]';
%!
%! % Randomization test to test whether the treatment sample comes from a
%! % population with mean of 56.2.
%! control = 56.2;
%! pval = randtest1 (treatment, control)
%! % The above is equivalent to:
%! % pval = randtest1 (treatment, control, 5000, @mean)


%!demo
%!
%! A = [21,26,33,22,18,25,26,24,21,25,35,28,32,36,38]';
%! GA = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8]';
%!
%! % Randomization test to test whether the sample A comes from a population
%! % population with mean of 30. Clusters of potentially correlated observations
%! % are defined in GA
%! M = 37;
%! pval = randtest1 ([A GA], M)
%! % The above is equivalent to:
%! % pval = randtest1 ([A GA], M, 5000, @mean)

%!test
%!
%! % Test various capabilities of randtest1
%! A = randn (3,1);
%! M = 0;
%! pval1 = randtest1 (A, M);
%! pval2 = randtest1 (A, M, []);
%! randtest1 (A, M, 500);
%! randtest1 (A, M, []);
%! A = randn (9,1);
%! pval3 = randtest1 (A, M, 5000);
%! pval4 = randtest1 (A, M, [], [], 1);
%! pval5 = randtest1 (A, M, [], 'mean', 1);
%! pval6 = randtest1 (A, M, [], @smoothmedian, 1);
