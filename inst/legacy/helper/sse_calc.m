function [SSb, SSw, K, g, MSb, MSw, dk] = sse_calc (x, groups, nvar)

  % Helper function file required for ibootci

  % Calculate error components of groups

  % Initialize
  gid = unique(groups);  % group ID
  K = numel(gid);        % number of groups
  n = numel(x{1});
  g = zeros(n,K);
  bSQ = zeros(K,nvar);
  wSQ = zeros(n,nvar);
  center = zeros(K,nvar);
  % Calculate within and between group variances
  for k = 1:K
    % Create group matrix
    g(:,k) = (groups == gid(k));
    for v = 1:nvar
      center(k,v) = sum(g(:,k) .* x{v}) / sum(g(:,k));
      wSQ(:,v) = wSQ(:,v) + g(:,k).*(x{v}-center(k,v)).^2;
    end
  end
  for v = 1:nvar
    bSQ(:,v) = (center(:,v) - mean(center(:,v))).^2;
  end
  SSb = sum(bSQ);         % Between-group SSE
  SSw = sum(wSQ);         % Within-group SSE
  g = logical(g);         % Logical array defining groups

  % Calculate mean squared error (MSE) and representative cluster size
  if nargout > 4
    nk = sum(g).';
    MSb = (sum(nk.*bSQ))/(K-1);
    MSw = SSw/(n-K);
    dk = mean(nk) - sum((sum(g)-mean(nk)).^2)/((K-1)*sum(g(:)));
  end

end
