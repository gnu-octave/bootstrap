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
%  algorithm uses a double bootstrap procedure to improve the accuracy of 
%  estimates for small-to-medium sample sizes [4-6].
%
%  stats = bootknife(data)
%  stats = bootknife(data,nboot)
%  stats = bootknife(data,nboot,bootfun)
%  stats = bootknife(data,nboot,bootfun,alpha)
%  stats = bootknife(data,nboot,bootfun,alpha,strata)
%  stats = bootknife(data,[2000,200],@mean,0.05,[])    % Default values
%  [stats,bootstat] = bootknife(...)
%  [stats,bootstat] = bootknife(...)
%  [stats,bootstat,bootsam] = bootknife(...)
%
%  stats = bootknife(data) resamples from the rows of a data sample (column 
%  vector or a matrix) and returns a column vector, which from top-to-
%  bottom contains the bootstrap bias-corrected estimate of the population 
%  mean [5,6], the bootstrap standard error of the mean, and calibrated 95% 
%  percentile bootstrap confidence intervals (lower and upper limits) [4]. 
%  Double bootstrap is used to improve the accuracy of the returned 
%  statistics, with the default number of outer (first) and inner (second) 
%  bootknife resamples being 2000 and 200 respectively. For confidence 
%  intervals, this is achieved by calibrating the lower and upper interval 
%  ends to have tail probabilities of 2.5% and 97.5% [4]. 
%
%  stats = bootknife(data,nboot) also specifies the number of bootknife 
%  samples. nboot can be a scalar, or vector of upto two positive integers. 
%  By default, nboot is [2000,200], which implements double bootstrap with 
%  the 2000 outer (first) and 200 inner (second) bootknife resamples. 
%  Larger number of resamples are recommended to reduce the Monte Carlo 
%  error, particularly of confidence intervals. If nboot provided is 
%  scalar, the value of nboot corresponds to the number of outer (first) 
%  bootknife resamples, and the default number of second bootstrap 
%  resamples (200) applies. Note that one can get away with a lower number 
%  of resamples in the second bootstrap (to reduce the computational 
%  expense of the double bootstrap) since the algorithm uses linear 
%  interpolation to achieve near asymptotic calibration of confidence 
%  intervals [3]. 
%    Setting the second element of nboot to 0 enforces a single bootstrap 
%  procedure. This can be useful if the purpose is just to get a quick, 
%  unbiased estimate of the bootstrap standard error by bootknife 
%  resampling. However, NaN will be returned for the confidence interval 
%  limits because uncalibrated percentile confidence intervals from a 
%  single bootstrap are considered too inaccurate.
%
%  stats = bootknife(data,nboot,bootfun) also specifies bootfun, a function 
%  handle (e.g. specified with @) or a string indicating the name of the 
%  function to apply to the data (and each bootknife resample). The default
%  value of bootfun is 'mean'.
%
%  stats = bootknife(data,nboot,bootfun,alpha) where alpha sets the lower 
%  and upper confidence interval ends to be 100 * (alpha/2)% and 100 * 
%  (1-alpha/2)% respectively. Central coverage of the intervals is thus 
%  100*(1-alpha)%. alpha should be a scalar value between 0 and 1. Default
%  is 0.05. If alpha is empty, NaN is returned for the confidence interval
%  ends.
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
%  statistics calculated over the (first or outer level of) bootknife 
%  resamples. 
%
%  [stats,bootstat,bootsam] = bootknife(...) also returns bootsam, the  
%  matrix of indices used for bootknife resampling. Each column in bootsam
%  corresponds to one bootknife sample and contains the row indices of the 
%  values drawn from the nonscalar data argument to create that sample.
%
%  Bibliography:
%  [1] Hesterberg T.C. (2004) Unbiasing the Bootstrap—Bootknife Sampling 
%        vs. Smoothing; Proceedings of the Section on Statistics & the 
%        Environment. Alexandria, VA: American Statistical Association.
%  [2] Davison et al. (1986) Efficient Bootstrap Simulation.
%        Biometrika, 73: 555-66
%  [3] Gleason, J.R. (1988) Algorithms for Balanced Bootstrap Simulations. 
%        The American Statistician. Vol. 42, No. 4 pp. 263-266
%  [4] Hall, Lee and Young (2000) Importance of interpolation when
%        constructing double-bootstrap confidence intervals. Journal
%        of the Royal Statistical Society. Series B. 62(3): 479-491
%  [5] Ouysee, R. (2011) Computationally efficient approximation for 
%        the double bootstrap mean bias correction. Economics Bulletin, 
%        AccessEcon, vol. 31(3), pages 2388-2403.
%  [6] Davison A.C. and Hinkley D.V (1997) Bootstrap Methods And Their 
%        Application. Chapter 3, pg. 104
%
%  bootknife v1.3.0.0 (16/05/2022)
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


function [stats, T1, idx] = bootknife (x, nboot, bootfun, alpha, strata)
  
  % Error checking
  if nargin < 1
    error ('data must be provided')
  end

  % Set defaults or check for errors
  if nargin < 2
    nboot = [2000, 200];
  else
    if isempty(nboot)
      nboot = [2000, 200];
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
  elseif ~isempty (alpha) 
    if size (strata, 1) ~= size (x, 1)
      error('strata must be a column vector with the same number of rows as the data')
    end
  end

  % Determine properties of the data (x)
  n = size (x, 1);

  % Initialize
  B = nboot(1);
  if (numel (nboot) > 1)
    C =  nboot(2);
  else
    C = 200;
  end
  T1 = zeros (1, B);
  idx = zeros (n, B);
  c = ones (n, 1) * B;

  % Perform balanced bootknife resampling
  if ~isempty (strata)
    % Get strata IDs
    gid = unique (strata);  % strata ID
    K = numel (gid);        % number of strata
    % Create strata matrix
    g = false (n,K);
    for k = 1:K
      g(:, k) = (strata == gid(k));
      [~, ~, idx(g(:, k), :)] = bootknife (x(g(:, k), :), [B, 0], bootfun, []);
      rows = find (g(:, k));
      idx(g(:, k), :) = rows(idx(g(:, k), :));
    end
  else
    for b = 1:B
      % Choose which rows of the data to sample
      r = b - fix ((b-1)/n) * n;
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
  for b = 1:B
    % Perform data sampling
    X = x(idx(:, b), :);
    % Function evaluation on bootknife sample
    T1(b) = feval (bootfun, X);
  end
 
  % Calculate the bootstrap standard error, bias and confidence intervals 
  % Bootstrap standard error estimation
  T0 = feval (bootfun,x);
  if C > 0
    U = zeros (1, B);
    M = zeros (1, B);
    V = zeros (1, B);
    % Iterated bootstrap resampling for greater accuracy
    for b = 1:B
      [~, T2] = bootknife (x(idx(:,b),:), [C,0], bootfun, [], strata);
      % Use quick interpolation to find the probability that T2 <= T0
      I = (T2 <= T0);
      u = sum (I);
      U(b) =  interp1q ([max([min(T2), max(T2(I))]);...
                         min([max(T2), min(T2(~I))])],...
                        [u; min(u+1, C)] / C,...
                        T0);
      if isnan (U(b))
        U(b) = u / C;
      end
      M(b) = mean (T2);
      V(b) = var (T2, 1);
    end
    % Double bootstrap bias estimation
    b = mean (T1) - T0;
    c = mean (M) - 2 * mean (T1) + T0;
    bias = b-c;
    % Double bootstrap standard error
    se = sqrt (var (T1, 1)^2 / mean (V));
    if ~isempty(alpha)
      % Calibrate tail probabilities to half of alpha
      l = quantile (U, [alpha/2, 1-alpha/2]);
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
    % (Uncalibrated) percentile confidence intervals are inaccurate and 
    % so we don't calculate them when a single bootstrap is requested
    ci = nan (1, 2);
  end
  
  % Prepare stats output argument
  stats = [T0 - bias; se; ci.'];
  
end
