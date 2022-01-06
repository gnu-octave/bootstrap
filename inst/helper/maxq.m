function q = maxq (Y,g,ref,bootfun,nvar)

  % Helper function file required for bootnhst

  % Calculate maximum difference between bootfun output of all the groups

  % Get size and of the data vector or matrix
  [m,n] = size(Y);
  N = size(g,1);
  if (nvar > 1)
    n = 1;
  end

  % Calculate the number (k) of unique groups
  gk = unique(g);
  k = numel(gk);

  % Calculate bootfun on the data from each group
  d = zeros(k,n);
  Var = zeros(k,1);
  nk = zeros(size(gk));
  for j = 1:k
    nk(j) = sum(g==gk(j));
    d(j,:) = feval(bootfun,Y(g==gk(j),:));
    Var(j,1) = ((nk(j)-1)/(N-k)) * jack(Y(g==gk(j),:), bootfun).^2;
  end
  nk_bar = sum(nk.^2)./sum(nk); % mean weighted sample size
  Var = sum(Var.*nk/nk_bar);    % pooled sampling variance weighted by sample size

  % Calculate weights to correct for unequal sample size  
  % when calculating standard error of the difference
  w = nk_bar./nk;

  % Calculate the q-ratio test statistic 
  if isempty(ref)
    % Calculate Tukey's q-ratio for maximum difference between  
    % bootfun for all sample pairwise comparisons
    %
    % Note that Tukey's q-ratio here does not have the sqrt(2) factor. 
    %
    % Bibliography:
    %  [1] https://en.wikipedia.org/wiki/Tukey%27s_range_test
    %  [2] https://cdn.graphpad.com/faq/1688/file/1688MulitpleComparisonAlgorithms(1).pdf
    %  [3] www.graphpad.com/guides/prism/latest/statistics/stat_the_methods_of_tukey_and_dunne.htm
    %
    [d,i] = sort(d,1);
    range = abs(d(k,:)-d(1,:));
    q = range / sqrt(Var * (w(i(k)) + w(i(1))));
  else
    % Calculate Dunnett's q-ratio for maximum difference between  
    % bootfun for test vs. control samples
    % Dunnett's q-ratio is similar to Student's t-statistic
    [range, i] = max(abs((d-ones(k,1)*d(ref,:))));
    q = range / sqrt(Var * (w(ref) + w(i)));
  end

end
