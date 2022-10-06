%  Function File: bootclust 
%
%  Two-stage nonparametric bootstrap sampling with shrinkage correction for
%  clustered data [1-4]. 
%
%  stats = bootclust (data,groups)
%  stats = bootclust (data,groups,nboot)
%  stats = bootclust (data,groups,nboot,alpha)
%  [stats, bootstat] = bootclust (...)
%  bootclust (data,groups,...)
%
%  stats = bootclust(data,groups) resamples both the group means of the data and 
%  their residuals. Bootstrap samples are formed by adding the resampled residuals
%  to the resampled means. Bootstrap statitistics are the mean of the bootstrap. 
%  The output structure, stats, contains the following fields:
%    original: contains the mean of the original data 
%    bias: contains the bootstrap estimate of bias
%    std_error: contains the bootstrap standard error
%    CI_lower: contains the lower bound of the bootstrap confidence interval
%    CI_upper: contains the upper bound of the bootstrap confidence interval
%  The confidence intervals are 95% percentile intervals.
%
%  stats = bootclust (data,groups,nboot) also specifies the number of bootstrap
%  samples. nboot must be a positive finite scalar. By default, nboot is 2000.
%
%  stats = bootclust (data,groups,nboot,alpha) where alpha sets the lower 
%  and upper bounds of the confidence interval(s). The value(s) in alpha must be 
%  between 0 and 1. If alpha is a scalar value, the nominal lower and upper
%  percentiles of the confidence are 100*(alpha/2)% and 100*(1-alpha/2)%
%  respectively, and nominal central coverage of the intervals is 100*(1-alpha)%.
%  If alpha is a vector with two elements, alpha becomes the quantiles for the
%  confidence intervals, and the intervals become percentile bootstrap confidence
%  intervals.
%
%  [stats, bootstat] = bootclust (...) also returns bootstat, a vector of
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
%
%  bootclust v0.5.0.0 (06/10/2022)
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
  
function stats = bootclust(data,groups,nboot,alpha)

  % Error checking
  if nargin < 2
    error('bootclust usage: ''bootci (data,groups)''; atleast 2 input arguments required');
  end
  if nargin < 3
    nboot = 2000;
  end
  if nargin < 4
    alpha = 0.05;
  end

  % Convert alpha to quantiles
  if numel(alpha) < 2
    % If the alpha provided is a scalar, convert it to quantiles
    q = [alpha / 2, 1 - alpha / 2];
  else
    q = alpha;
  end

  % Calculate the means and residuals for each group/cluster
  [mu,resid,K,g] = clustmean(data,groups);

  % Redefine bootfun for two-stage balanced, bootknife resampling
  bootfun = @(resid) clustboot(mu,resid,K,g);

  % Run resampling and calculation of bootstrao statistics
  if nargout < 1
    bootknife(resid, nboot, bootfun, q);
  else
    [stats, bootstat] = bootknife(resid, nboot, bootfun, q);
  end

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
