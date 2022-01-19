function q = maxq (Y,g,bootfun,ref)

  % Helper function file required for bootnhst

  % Calculate maximum studentized difference between bootfun output of all the groups

  % Get size and of the data vector or matrix
  [m,nvar] = size(Y);
  N = size(g,1);

  % Calculate the number (k) of unique groups
  gk = unique(g);
  k = numel(gk);

  % Compute the estimate (theta) and it's pooled (weighted mean) sampling variance 
  theta = zeros(k,1);
  SE = zeros(k,1);
  Var = zeros(k,1);
  B = 200;
  t = zeros(B,1); 
  nk = zeros(size(gk));
  for j = 1:k
    nk(j) = sum(g==gk(j));
    theta(j) = feval(bootfun,Y(g==gk(j),:));
    % Compute unbiased estimate of the standard error by bootknife resampling
    if nvar > 1
      t = zeros(B,1); 
      for b = 1:B
        idx = 1+fix(rand(nk(j)-1,1)*nk(j));
        tmp = Y(g==gk(j),:);
        t(b) = feval(bootfun,tmp(idx,:));
      end
    else
      % Vectorized if data is univariate
      idx = 1+fix(rand(nk(j)-1,B)*nk(j));
      tmp = Y(g==gk(j),:);
      t = feval(bootfun,tmp(idx));
    end
    SE(j) = std(t);
    Var(j) = ((nk(j)-1)/(N-k)) * SE(j)^2;
  end
  nk_bar = sum(nk.^2)./sum(nk);  % weighted mean sample size
  Var = sum(Var.*nk/nk_bar);     % pooled sampling variance weighted by sample size

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
    %  [2] https://cdn.graphpad.com/faq/1688/file/MulitpleComparisonAlgorithmsPrism8.pdf
    %  [3] www.graphpad.com/guides/prism/latest/statistics/stat_the_methods_of_tukey_and_dunne.htm
    %
    [theta,i] = sort(theta,1);
    range = abs(theta(k) - theta(1));
    q = range / sqrt(Var * (w(i(k)) + w(i(1))));
  else
    % Calculate Dunnett's q-ratio for maximum difference between  
    % bootfun for test vs. control samples
    % Dunnett's q-ratio is similar to Student's t-statistic
    [range, i] = max(abs((theta - ones(k,1) * theta(ref))));
    q = range / sqrt(Var * (w(ref) + w(i)));
  end
  
end
