%  Function File: bootknife
%
%  Bootknife resampling  
%
%  This function takes a data sample (containing n rows) and uses bootstrap 
%  techniques to calculate a bias-corrected parameter estimate, a standard 
%  error, and 95% confidence intervals. Specifically, the method uses 
%  bootknife resampling [1], which involves creating leave-one-out 
%  jackknife samples of size n - 1 and then drawing samples of size n with 
%  replacement from the jackknife samples. The resampling of data rows is 
%  balanced in order to reduce Monte Carlo error [2,3]. By default, the 
%  bootstrap confidence intervals are bias-corrected and accelerated (BCa) 
%  [4-5]. BCa intervals are fast to compute and have good coverage and 
%  correctness when combined with bootknife resampling as it is here [1], 
%  but it may not have the intended coverage when sample size gets very 
%  small. If double bootstrap is requested, the algorithm uses calibration 
%  to improve the accuracy of bias-corrected estimates and confidence 
%  intervals for small-to-medium sample sizes [6-8]. 
%
%  stats = bootknife(data)
%  stats = bootknife(data,nboot)
%  stats = bootknife(data,nboot,bootfun)
%  stats = bootknife(data,nboot,{bootfun,bootfun_args})
%  stats = bootknife(data,nboot,bootfun,alpha)
%  stats = bootknife(data,nboot,bootfun,alpha,strata)
%  stats = bootknife(data,[2000,0],@mean,0.05,[])      % Default values
%  [stats,bootstat] = bootknife(...)
%  [stats,bootstat] = bootknife(...)
%  [stats,bootstat,bootsam] = bootknife(...)
%
%  stats = bootknife(data) resamples from the rows of a data sample (column 
%  vector or a matrix) and returns a column vector or matrix, whose rows  
%  (from top-to-bottom) contain the bootstrap bias-corrected estimate of   
%  the population mean(s) [7-8], the bootknife standard error of the mean  
%  [1], and 95% bias-corrected and accelerated (BCa) confidence intervals 
%  [1,4-5,9]. The data cannot contain values that are NaN or infinite.
%
%  stats = bootknife(data,nboot) also specifies the number of bootknife 
%  samples. nboot can be a scalar, or vector of upto two positive integers. 
%  By default, nboot is [2000,0], which implements a single bootstrap with 
%  the 2000 resamples, but larger numbers of resamples are recommended to  
%  reduce the Monte Carlo error, particularly for confidence intervals. If  
%  the second element of nboot is > 0, then the first and second elements  
%  of nboot correspond to the number of outer (first) and inner (second) 
%  bootknife resamples respectively. Double bootstrap is used to improve 
%  the accuracy of the bias-corrected estimate(s) and the confidence 
%  intervals. For confidence intervals, this is achieved by calibrating 
%  the lower and upper interval ends to have tail probabilities of 2.5% 
%  and 97.5% [5]. Note that one can get away with a lower number of 
%  resamples in the second bootstrap to reduce the computational expense 
%  of the double bootstrap (e.g. [2000,200]), since the algorithm uses 
%  linear interpolation to achieve near-asymptotic calibration of 
%  confidence intervals [3]. The confidence intervals calculated (with 
%  either single or double bootstrap) are transformation invariant and 
%  have more accuracy and correctness compared to intervals derived from 
%  normal theory or to simple percentile bootstrap confidence intervals. 
%
%  stats = bootknife(data,nboot,bootfun) also specifies bootfun, a function 
%  handle (e.g. specified with @), a string indicating the name of the 
%  function to apply to the data (and each bootknife resample), or a cell 
%  array where the first cell is the function handle or string, and other 
%  cells being arguments for that function, where the function must take 
%  data for the first input argument. bootfun can return a scalar value or 
%  vector. The default value(s) of bootfun is/are the (column) mean(s).
%  When bootfun is @mean or 'mean', residual narrowness bias of central 
%  coverage is eliminated by using the Student's t-distribution to expand 
%  the percentiles before applying the BCa adjustments as described in [9].
%    Note that bootfun MUST calculate a statistic representative of the 
%  finite data sample, it should NOT be an estimate of a population 
%  parameter. For example, for the variance, set bootfun to {@var,1}, not 
%  @var or {@var,0}. Smooth functions of the data are preferable, (e.g. use 
%  smoothmedian function instead of ordinary median). 
%
%  stats = bootknife(data,nboot,bootfun,alpha) where alpha sets the lower 
%  and upper confidence interval ends to be 100 * (alpha/2)% and 100 * 
%  (1-alpha/2)% respectively. Central coverage of the intervals is thus 
%  100*(1-alpha)%. alpha should be a scalar value between 0 and 1. If alpha 
%  is empty, NaN is returned for the confidence interval ends. The default
%  alpha is 0.05. 
%
%  stats = bootknife(data,nboot,bootfun,alpha,strata) also sets strata, 
%  which are numeric identifiers that define the grouping of the data rows
%  for stratified bootknife resampling. strata should be a column vector 
%  the same number of rows as the data. When resampling is stratified, 
%  the groups (or stata) of data are equally represented across the 
%  bootknife resamples. If this input argument is not specified or is 
%  empty, no stratification of resampling is performed. 
%
%  [stats,bootstat] = bootknife(...) also returns bootstat, a vector of
%  statistics calculated over the (first, or outer level of) bootknife 
%  resamples. 
%
%  [stats,bootstat,bootsam] = bootknife(...) also returns bootsam, the  
%  matrix of indices used for the (first, or outer level of) bootknife 
%  resampling. Each column in bootsam corresponds to one bootknife 
%  resample and contains the row indices of the values drawn from the 
%  nonscalar data argument to create that sample.
%
%  Bibliography:
%  [1] Hesterberg T.C. (2004) Unbiasing the Bootstrap—Bootknife Sampling 
%        vs. Smoothing; Proceedings of the Section on Statistics & the 
%        Environment. Alexandria, VA: American Statistical Association.
%  [2] Davison et al. (1986) Efficient Bootstrap Simulation.
%        Biometrika, 73: 555-66
%  [3] Gleason, J.R. (1988) Algorithms for Balanced Bootstrap Simulations. 
%        The American Statistician. Vol. 42, No. 4 pp. 263-266
%  [4] Efron (1987) Better Bootstrap Confidence Intervals. JASA, 
%        82(397): 171-185 
%  [5] Efron, and Tibshirani (1993) An Introduction to the
%        Bootstrap. New York, NY: Chapman & Hall
%  [6] Hall, Lee and Young (2000) Importance of interpolation when
%        constructing double-bootstrap confidence intervals. Journal
%        of the Royal Statistical Society. Series B. 62(3): 479-491
%  [7] Ouysee, R. (2011) Computationally efficient approximation for 
%        the double bootstrap mean bias correction. Economics Bulletin, 
%        AccessEcon, vol. 31(3), pages 2388-2403.
%  [8] Davison A.C. and Hinkley D.V (1997) Bootstrap Methods And Their 
%        Application. Chapter 3, pg. 104
%  [9] Hesterberg, Tim (2014), What Teachers Should Know about the 
%        Bootstrap: Resampling in the Undergraduate Statistics Curriculum, 
%        http://arxiv.org/abs/1411.5279
%
%  bootknife v1.4.5.0 (26/05/2022)
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


function [stats, T1, idx] = bootknife (x, nboot, bootfun, alpha, strata, idx)
  
  % Error checking
  if nargin < 1
    error ('data must be provided')
  end

  % Set defaults or check for errors
  if nargin < 2
    nboot = [2000, 0];
  else
    if isempty(nboot)
      nboot = [2000, 0];
    else
      if ~isa (nboot, 'numeric')
        error('nboot must be numeric');
      end
      if any (nboot ~= abs (fix (nboot)))
        error ('nboot must contain positive integers')
      end    
      if numel (nboot)>3
        error ('nboot cannot contain more than 2 values')
      end
    end
  end
  if nargin < 3
    bootfun = @mean;
  else
    if iscell(bootfun)
      func = bootfun{1};
      args = bootfun(2:end);
      bootfun = @(x) feval(func, x, args{:});
    end
    if ischar (bootfun)
      % Convert character string of a function name to a function handle
      bootfun = str2func (bootfun);
    end
    if ~isa (bootfun, 'function_handle')
      error('bootfun must be a function name or function handle');
    end
  end
  if nargin < 4
    alpha = 0.05;
  elseif ~isempty (alpha) 
    if ~isa (alpha,'numeric') || numel (alpha)~=1
      error('alpha must be a numeric scalar value');
    end
    if (alpha <= 0) || (alpha >= 1)
      error('alpha must be a value between 0 and 1');
    end
  end
  if nargin < 5
    strata = [];
  elseif ~isempty (strata) 
    if size (strata, 1) ~= size (x, 1)
      error('strata must be a column vector with the same number of rows as the data')
    end
  end

  % Determine properties of the data (x)
  n = size (x, 1);

  % Set number of outer and inner bootknife resamples
  B = nboot(1);
  if (numel (nboot) > 1)
    C =  nboot(2);
  else
    C = 0;
  end

  % Evaluate bootfun on the data
  T0 = feval (bootfun, x);
  if all (size (T0) > 1)
    error ('bootfun must return either a scalar or a vector')
  end

  % If the function of the data is a vector, calculate the statistics for each element 
  m = numel(T0);
  if m > 1
    % Evaluate bootknife for each element of the output of bootfun
    stats = zeros (4, m);
    T1 = zeros (m, B);
    for j = 1:m
      out = @(x, j) x(j);
      func = @(x) out (bootfun(x), j); 
      if j > 1
        [stats(:,j), T1(j,:)] = bootknife (x, nboot, func, alpha, strata, idx);
      else
        [stats(:,j), T1(j,:), idx] = bootknife (x, nboot, func, alpha, strata);
      end
    end
    return
  end

  % Evaluate strata input argument
  if ~isempty (strata)
    % Get strata IDs
    gid = unique (strata);  % strata ID
    K = numel (gid);        % number of strata
    % Create strata matrix
    g = false (n,K);
    for k = 1:K
      g(:, k) = (strata == gid(k));
    end
  end

  % Perform balanced bootknife resampling
  if nargin < 6
    % Initialize
    T1 = zeros (1, B);
    idx = zeros (n, B);
    c = ones (n, 1) * B;
    % Calculate row indices for resampling
    if ~isempty (strata)
      for k = 1:K
        [~, ~, idx(g(:, k), :)] = bootknife (x(g(:, k), :), [B, 0], bootfun, []);
        rows = find (g(:, k));
        idx(g(:, k), :) = rows(idx(g(:, k), :));
      end
    else
      for b = 1:B
        % Choose which rows of the data to sample
        r = b - fix ((b - 1) / n) * n;
        for i = 1:n
          d = c;   
          d(r) = 0;
          if ~sum (d)
            d = c;
          end
          j = sum ((rand (1) >= cumsum (d ./sum (d)))) + 1;
          idx(i, b) = j;
          c(j) = c(j) - 1;
        end 
      end
    end
  end
  for b = 1:B
    % Perform data sampling
    X = x(idx(:, b), :);
    % Function evaluation on bootknife sample
    T1(b) = feval (bootfun, X);
  end
 
  % Calculate the bootstrap bias, standard error and confidence intervals 
  if C > 0
    U = zeros (1, B);
    M = zeros (1, B);
    V = zeros (1, B);
    % Iterated bootstrap resampling for greater accuracy
    for b = 1:B
      [~, T2] = bootknife (x(idx(:, b), :), [C, 0], bootfun, [], strata);
      % Use quick interpolation to find the probability that T2 <= T0
      I = (T2 <= T0);
      u = sum (I);
      t2 = [max([min(T2), max(T2(I))]),...
            min([max(T2), min(T2(~I))])];
      % Linear interpolation
      U(b) = ((t2(2) - T0) * u / C + (T0 - t2(1)) * min ((u + 1) / C, 1)) /...
            (t2(2) - t2(1));
      if isnan (U(b))
        U(b) = u / C;
      end
      M(b) = mean (T2);
      V(b) = var (T2, 1);
    end
    % Double bootstrap bias estimation
    b = mean (T1) - T0;
    c = mean (M) - 2 * mean (T1) + T0;
    bias = b - c;
    % Bootstrap standard error
    se = std (T1, 1);
    if ~isempty(alpha)
      % Calibrate tail probabilities to half of alpha
      l = quantile (U, [alpha / 2, 1 - alpha / 2]);
      % Calibrated percentile bootstrap confidence intervals
      ci = quantile (T1, l);
    else
      ci = nan (1, 2);
    end
  else
    % Bootstrap bias estimation
    bias = mean (T1) - T0;
    % Bootstrap standard error
    se = std (T1, 1);
    if ~isempty(alpha)
      % Bias-corrected and accelerated bootstrap confidence intervals 
      % Create distribution functions
      stdnormcdf = @(x) 0.5 * (1 + erf (x / sqrt (2)));
      stdnorminv = @(p) sqrt (2) * erfinv (2 * p-1);
      % Calculate the median bias correction z0
      z0 = stdnorminv (sum (T1 < T0)/ B);
      if isinf (z0) || isnan (z0)
        error('unable to calculate the bias correction z0')
      end
      % Use the Jackknife to calculate the acceleration constant
      T = zeros (n,1);
      for i = 1:n
        T(i) = feval (bootfun, x(1:end ~= i, :));
      end
      % Calculate empirical influence function
      if ~isempty(strata)
        gk = sum (g .* repmat (sum (g), n, 1), 2);
        U = (gk - 1) .* (mean (T) - T);   
      else
        U = (n - 1) * (mean (T) - T);     
      end
      a = sum (U.^3) / (6 * sum (U.^2) ^ 1.5);
      if strcmp (func2str (bootfun), 'mean')
        % If bootfun is the mean, expand percentiles using Student's 
        % t-distribution to improve central coverage for small samples
        studinv = @(p, df) - sqrt ( df ./ betaincinv (2 * p, df / 2, 0.5) - df);
        alpha = stdnormcdf (studinv (alpha / 2, n - 1)) * 2;      
      end
      % Calculate BCa percentiles
      z1 = stdnorminv(alpha / 2);
      z2 = stdnorminv(1 - alpha / 2);
      l = cat(2, stdnormcdf (z0 + ((z0 + z1) / (1 - a * (z0 + z1)))),... 
                 stdnormcdf (z0 + ((z0 + z2) / (1 - a * (z0 + z2)))));
      ci = quantile (T1, l);
    else
      ci = nan (1, 2);
    end
  end
  
  % Prepare stats output argument
  stats = [T0 - bias; se; ci.'];
  
end

