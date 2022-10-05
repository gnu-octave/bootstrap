% function ci = bootclust (data,groups,nboot,alpha)
%
% Two-stage nonparametric bootstrap sampling with shrinkage correction for
% clustered data [1-4]. groups should specify a column vector (or matrix) of
% numeric identifiers with the same number of rows as the data. 
%
% By resampling residuals, this bootstrap method can be used when cluster sizes
% are unequal. However, cluster samples are assumed to be taken from populations
% with equal variance.
%
% References:
%  [1] Davison and Hinkley (1997) Bootstrap Methods and their
%       application. Chapter 3: pg 97-100
%  [2] Ng, Grieve and Carpenter (2013) The Stata Journal.
%       13(1): 141-164
%  [3] Gomes et al. (2012) Medical Decision Making. 32(2): 350-361
%  [4] Gomes et al. (2012) Health Econ. 21(9):1101-18
  
function stats = bootclust(data,groups,nboot,alpha)

  % Error checking
  if nargin < 3
    nboot = 2000;
  end
  if nargin < 4
    alpha = 0.05;
  end

  % Convert alpha to quantiles
  q = [alpha / 2, 1 - alpha / 2];

  % Calculate the means and residuals for each group/cluster
  [mu,resid,K,g] = clustmean(data,groups);

  % Redefine bootfun for two-stage balanced, bootknife resampling
  bootfun = @(resid) clustboot(mu,resid,K,g);

  % Run resampling and calculation of bootstrao statistics
  if nargout < 1
    bootknife(resid, nboot, bootfun, q);
  else
    stats = bootknife(resid, nboot, bootfun, q);
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
