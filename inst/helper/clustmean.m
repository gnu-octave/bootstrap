function [mu, Z, K, g] = clustmean (x, clusters, nvar)

  % Helper function file required for ibootci

  % Calculates shrunken cluster means and residuals for cluster bootstrap
  % See also bootclust function below

  % Calculate sum-of-squared error components and number of clusters
  [SSb, SSw, K, g] = sse_calc (x, clusters, nvar);
  SSb = sum(SSb);
  SSw = sum(SSw);

  % Calculate cluster means in the original sample
  mu = cell(1,nvar);
  for v = 1:nvar
    for k = 1:K
      mu{v}(k,:) = mean(x{v}(g(:,k),:));
    end
  end

  % Calculate shrunken cluster means from the original sample
  nk = sum(g).';
  dk = mean(nk) - sum((sum(g)-mean(nk)).^2)/((K-1)*sum(g(:)));
  c = 1 - sqrt(max(0,(K/(K-1)) - (SSw./(dk.*(dk-1).*SSb))));
  for v = 1:nvar
    mu{v} = bsxfun(@plus, c*mean(mu{v}),(1-c)*mu{v});
  end

  % Calculate residuals from the sample and cluster means
  Z = cell(1,nvar);
  for v = 1:nvar
    for k = 1:K
      Z{v}(g(:,k),:) = bsxfun(@minus, x{v}(g(:,k),:), mu{v}(k,:));
      Z{v}(g(:,k),:) = Z{v}(g(:,k),:) ./ sqrt(1-dk^-1);
    end
  end

end
