%  Function File: bootclust 
%
%  Two-stage nonparametric bootstrap sampling with shrinkage correction for
%  the mean of clustered data [1-4]. 
%
%  ci = bootclust (data,groups)
%  ci = bootclust (data,groups,nboot)
%  ci = bootclust (data,groups,nboot,alpha)
%  [ci, bootstat] = bootclust (...)
%  bootclust (data,groups,...)
%
%  ci = bootclust(data,groups) resamples both the group means of the data and 
%  their residuals. Bootstrap samples are formed by adding the resampled residuals
%  to the resampled means. Bootstrap statitistics are the means of the bootstrap 
%  samples. The confidence intervals returned are expanded percentile intervals [5].
%
%  ci = bootclust (data,groups,nboot) also specifies the number of bootstrap
%  samples. nboot must be a positive finite scalar. By default, nboot is 2000.
%
%  ci = bootclust (data,groups,nboot,alpha) where alpha is a scalar value that
%  sets the lower and upper bounds of the confidence interval(s). The value of
%  alpha must be between 0 and 1 and the nominal central coverage of the intervals
%  is 100*(1-alpha)%. The default value for alpha is 0.05. 
%
%  [ci, bootstat] = bootclust (...) also returns bootstat, a vector of
%  statistics calculated over the bootstrap samples.
%
%  bootclust(data,...); returns a pretty table of the output including
%  the bootstrap settings and the result of evaluating bootfun on the
%  data along with bootstrap estimates of bias, standard error, and
%  lower and upper 100*(1-alpha)% confidence limits.
%
% References:
%  [1] Davison and Hinkley (1997) Bootstrap Methods and their
%       application. Chapter 3: pg 97-100
%  [2] Ng, Grieve and Carpenter (2013) The Stata Journal.
%       13(1): 141-164
%  [3] Gomes et al. (2012) Medical Decision Making. 32(2): 350-361
%  [4] Gomes et al. (2012) Health Econ. 21(9):1101-18
%  [5] Hesterberg, Tim (2014), What Teachers Should Know about the 
%        Bootstrap: Resampling in the Undergraduate Statistics Curriculum, 
%        http://arxiv.org/abs/1411.5279
%
%  bootclust (version 2022.10.08)
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
  
function [ci, bootstat] = bootclust(data,groups,nboot,alpha)

  % Error checking
  if nargin < 2
    error('bootclust usage: ''bootclust (data,groups)''; atleast 2 input arguments required');
  end
  if size (groups, 1) ~= size (data, 1)
    error ('bootclust: group should be a column vector or cell array with the same number of rows as the data');
  end
  if any(isnan(data)) || any(isinf(data))
    error ('bootclust: data must contain finite numeric values');
  end
  if (nargin < 3)
    nboot = 2000;
  else
    if isempty(nboot)
      nboot = 2000;
    else
      if ~isa (nboot, 'numeric')
        error ('bootclust: nboot must be numeric');
      end
      if any (nboot ~= abs (fix (nboot)))
        error ('bootclust: nboot must contain positive integers');
      end    
      if (numel (nboot) > 1)
        error ('bootclust: nboot must be scalar');
      end
    end
  end
  if (nargin < 4)
    alpha = 0.05;
  elseif ~isempty (alpha) 
    if ~isa (alpha,'numeric') 
      error ('bootclust: alpha must be numeric');
    end
    if any ((alpha < 0) | (alpha > 1))
      error ('bootclust: value(s) in alpha must be between 0 and 1');
    end
    if numel(alpha) > 1
      error ('bootclust: alpha must be a scalar value')
    end
  end
  % Evaluate group input argument
  if ~isnumeric (groups)
    % Convert group to numeric ID
    [junk1,junk2,groups] = unique (groups);
    clear junk1 junk2;
  end

  % Calculate the shrunken means and residuals for each group/cluster
  [mu,resid,K,g] = clustmean (data,groups);

  % Calculate expanded percentiles
  stdnormcdf = @(x) 0.5 * (1 + erf (x / sqrt (2)));
  stdnorminv = @(p) sqrt (2) * erfinv (2 * p-1);
  if exist('betaincinv','file')
    studinv = @(p, df) - sqrt ( df ./ betaincinv (2 * p, df / 2, 0.5) - df);
  else
    % Earlier versions of matlab do not have betaincinv
    % Instead, use betainv from the Statistics and Machine Learning Toolbox
    try 
      studinv = @(p, df) - sqrt ( df ./ betainv (2 * p, df / 2, 0.5) - df);
    catch
      % Use the Normal distribution (i.e. do not expand quantiles) if
      % either betaincinv or betainv are not available
      studinv = @(p,df) sqrt (2) * erfinv (2 * p-1);
      warning ('bootclust could not expand the percentiles; interval coverage will be poor when thenumber of clusters is small')
    end
  end
  adj_alpha = stdnormcdf (studinv (alpha / 2, K - 1)) * 2;
  q = [adj_alpha / 2, 1 - adj_alpha / 2];

  % Redefine bootfun for two-stage balanced, bootknife resampling
  bootfun = @(resid) clustboot (mu,resid,K,g);

  % Run resampling and calculation of bootstrao statistics
  bootstat = bootfun (boot (resid, nboot, false));

  % Calculate percentile confidence intervals
  [cdf, t1] = empcdf (bootstat, 1);
  ci = arrayfun ( @(p) interp1 (cdf, t1, p, 'linear'), q);

  % Transpose result
  ci = ci';
  bootstat = bootstat';

end

%--------------------------------------------------------------------------

function [mu, resid, K, g] = clustmean (x, clusters)

  % Calculate shrunken cluster means and residuals for cluster bootstrap resampling

  % Calculate sum-of-squared error components and number of clusters
  [SSb, SSw, K, g] = sse_calc (x, clusters);
  SSb = sum(SSb);
  SSw = sum(SSw);

  % Calculate cluster means in the original sample
  mu = zeros(K,1);
  for k = 1:K
    mu(k) = mean(x(g(:,k),:));
  end

  % Calculate shrunken cluster means from the original sample
  nk = sum(g).';
  dk = mean(nk) - sum((sum(g)-mean(nk)).^2)/((K-1)*sum(g(:)));
  c = 1 - sqrt(max(0,(K/(K-1)) - (SSw./(dk.*(dk-1).*SSb))));
  mu = bsxfun(@plus, c*mean(mu),(1-c)*mu);

  % Calculate residuals from the sample and cluster means
  resid = zeros(sum(nk),1);
  for k = 1:K
    resid(g(:,k),:) = bsxfun(@minus, x(g(:,k),:), mu(k,:));
    resid(g(:,k),:) = resid(g(:,k),:) ./ sqrt(1-dk^-1);
  end

end

%--------------------------------------------------------------------------

function [SSb, SSw, K, g] = sse_calc (x, groups)

  % Calculate error components of groups

  % Initialize
  gid = unique(groups);  % group ID
  K = numel(gid);        % number of groups
  n = numel(x);
  g = zeros(n,K);
  bSQ = zeros(K,1);
  wSQ = zeros(n,1);
  center = zeros(K,1);
  % Calculate within and between group variances
  for k = 1:K
    % Create group matrix
    g(:,k) = (groups == gid(k));
    center(k,1) = sum(g(:,k) .* x) / sum(g(:,k));
    wSQ(:,1) = wSQ(:,1) + g(:,k).*(x-center(k,1)).^2;
  end
  bSQ(:,1) = (center(:,1) - mean(center(:,1))).^2;
  SSb = sum(bSQ);         % Between-group SSE
  SSw = sum(wSQ);         % Within-group SSE
  g = logical(g);         % Logical array defining groups

end

%--------------------------------------------------------------------------

function T = clustboot (mu, resid, K, g)

  % The cluster means are resampled and combined with the residuals before
  % function evaluation

  % Calculate data dimensions
  [n,nboot] = size(resid);

  % Balanced, bootknife resampling of cluster means
  bootmu = boot(mu,nboot,true);

  % Combine residuals with resampled cluster means
  X = zeros(n,nboot);
  for k = 1:K
    X(g(:,k),:) = bsxfun(@plus, resid(g(:,k),:), bootmu(k,:));
  end

  % Calculate bootstrap statistic(s)
  T = mean(X);

end

%--------------------------------------------------------------------------

function [F, x] = empcdf (bootstat, c)

  % Subfunction to calculate empirical cumulative distribution function of bootstat
  %
  % Set c to:
  %  1 to have a complete distribution with F ranging from 0 to 1
  %  0 to avoid duplicate values in x
  %
  % Unlike ecdf, empcdf uses a denominator of N+1

  % Check input argument
  if ~isa(bootstat,'numeric')
    error ('bootknife:empcdf: bootstat must be numeric');
  end
  if all(size(bootstat)>1)
    error ('bootknife:empcdf: bootstat must be a vector');
  end
  if size(bootstat,2)>1
    bootstat = bootstat.';
  end

  % Create empirical CDF
  bootstat = sort(bootstat);
  N = sum(~isnan(bootstat));
  [x,F] = unique(bootstat,'rows','last','legacy');
  F = F/(N+1);

  % Apply option to complete the CDF
  if c > 0
    x = [x(1);x;x(end)];
    F = [0;F;1];
  end

  % Remove impossible values
  F(isnan(x)) = [];
  x(isnan(x)) = [];
  F(isinf(x)) = [];
  x(isinf(x)) = [];

end

%--------------------------------------------------------------------------

%!test
%! data = [0.125, 0.127, 0.125, 0.126, 0.128, 0.118, 0.122, 0.12, 0.124, ...
%!         0.119, 0.123, 0.125, 0.125, 0.124, 0.126, 0.126, 0.128, 0.126, ...
%!         0.127, 0.129, 0.118, 0.129, 0.127, 0.12, 0.121]';
%! groups = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, ...
%!           4, 4, 4, 4, 4, 5, 5, 5, 5, 5]';
%! boot (1, 1, false, [], 1); # Set random seed
%! ci = bootclust (data,groups);
