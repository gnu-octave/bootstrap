function T = maxdiff(Y,g,ref,bootfun,nvar)

  % Helper function file required for bootnhst

  % Calculate maximum difference between bootfun output of all the groups

  % Get size and of the data vector or matrix
  [m,n] = size(Y);
  if (nvar > 1)
    n = 1;
  end

  % Calculate the number (k) of unique groups
  gk = unique(g);
  k = numel(gk);

  % Calculate bootfun on the data from each group
  Z = zeros(k,n);
  for j = 1:k
    Z(j,:) = feval(bootfun,Y(g==gk(j),:));
  end

  % Calculate maximum difference test statistic (T) 
  % This statistic is simpler and more intuitive than calculating F 
  % or MSE, and makes it more straightforward to calculate post-tests
  if isempty(ref)
    % Maximum difference for all pairwise comparisons
    Z = sort(Z,1);
    T = Z(k,:)-Z(1,:); % sign always positive
  else
    % Maximum absolute difference vs. reference/control
    T = max(abs(Z-ones(k,1)*Z(ref,:)));
  end

end
