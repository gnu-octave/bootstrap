function maxT = maxstat (Y, g, nboot, bootfun, ref, clusters, strata)

  % Helper function file required for ibootnhst

  % Calculate maximum test statistic

  % Get size and of the data vector or matrix
  [m,nvar] = size(Y);

  % Get data structure information
  if isempty(clusters)
    N = size(g,1);
  else
    N = numel(unique(clusters));
  end
  if isempty(strata)
    l = 1;
  else
    l = numel(unique(strata)); % number of strata
  end

  % Calculate the number (k) of unique groups
  gk = unique(g);
  k = numel(gk);

  % Compute the estimate (theta) and it's pooled (weighted mean) sampling variance 
  theta = zeros(k,1);
  SE = zeros(k,1);
  Var = zeros(k,1);
  t = zeros(nboot,1); 
  nk = zeros(size(gk));
  for j = 1:k
    theta(j) = feval(bootfun,Y(g==gk(j),:));
    if ~isempty(clusters)
      % Compute unbiased estimate of the standard error by cluster-jackknife resampling
      opt = struct;
      opt.clusters = clusters(g==gk(j));
      nk(j) = numel(unique(opt.clusters));
      SE(j) = jack (Y(g==gk(j),:), bootfun, [], opt);
    elseif (nboot == 0)
      % If requested, compute unbiased estimates of the standard error using jackknife resampling
      nk(j) = sum(g==gk(j));
      SE(j) = jack (Y(g==gk(j),:), bootfun);
    else
      % Compute unbiased estimate of the standard error by bootknife resampling
      % Bootknife resampling involves less computation than Jackknife when sample sizes get larger
      nk(j) = sum(g==gk(j));
      if nvar > 1
        t = zeros(nboot,1); 
        for b = 1:nboot
          idx = 1+fix(rand(nk(j)-1,1)*nk(j));
          tmp = Y(g==gk(j),:);
          t(b) = feval(bootfun,tmp(idx,:));
        end
      else
        % Vectorized if data is univariate
        idx = 1+fix(rand(nk(j)-1,nboot)*nk(j));
        tmp = Y(g==gk(j),:);
        t = feval(bootfun,tmp(idx));
      end
      SE(j) = std(t);
    end
    Var(j) = ((nk(j)-1)/(N-k-(l-1))) * SE(j)^2;
  end
  if any(nk <= 1)
    error('the number of observations or clusters per group must be greater than 1')
  end
  nk_bar = sum(nk.^2)./sum(nk);  % weighted mean sample size
  Var = sum(Var.*nk/nk_bar);     % pooled sampling variance weighted by sample size

  % Calculate weights to correct for unequal sample size  
  % when calculating standard error of the difference
  w = nk_bar./nk;

  % Calculate the maximum test statistic 
  if isempty(ref)
    % Calculate Tukey-Kramer test statistic (without sqrt(2) factor)
    %
    % Bibliography:
    %  [1] https://en.wikipedia.org/wiki/Tukey%27s_range_test
    %  [2] https://cdn.graphpad.com/faq/1688/file/MulitpleComparisonAlgorithmsPrism8.pdf
    %  [3] www.graphpad.com/guides/prism/latest/statistics/stat_the_methods_of_tukey_and_dunne.htm
    idx = logical(triu(ones(k,k),1));
    i = (1:k)' * ones(1,k);
    j = ones(k,1) * (1:k);
    t = abs(theta(i(idx)) - theta(j(idx))) ./ sqrt(Var * (w(i(idx)) + w(j(idx))));;
  else
    % Calculate Dunnett's test statistic 
    t = abs((theta - theta(ref))) ./ sqrt(Var * (w + w(ref)));
  end
  maxT = max(t);
  
end
