function maxT = maxstat (Y, g, nboot, bootfun, ref, ISOCTAVE)

  % Helper function file required for bootnhst
  % Calculate maximum test statistic
  
  % maxstat cannot be a subfunction or nested function since 
  % Octave parallel threads won't be able to find it

  % Calculate the size of the data (N) and the number (k) of unique groups
  N = size(g,1);
  gk = unique(g);
  k = numel(gk);

  % Compute the estimate (theta) and it's pooled (weighted mean) sampling variance 
  theta = zeros(k,1);
  SE = zeros(k,1);
  Var = zeros(k,1);
  nk = zeros(size(gk));
  for j = 1:k
    if (nboot == 0)
      nk(j) = sum(g==gk(j));
      if strcmp (func2str(bootfun), 'mean')
        theta(j) = mean(Y(g==gk(j),:));
        % Quick calculation for the standard error of the mean
        SE(j) = std(Y(g==gk(j),:),0) / sqrt(nk(j));
      else
        theta(j) = bootfun(Y(g==gk(j),:));
        % If requested, compute unbiased estimates of the standard error using jackknife resampling
        SE(j) = jack(Y(g==gk(j),:), bootfun);
      end
    else
      % Compute unbiased estimate of the standard error by balanced bootknife resampling
      % Bootknife resampling involves less computation than Jackknife when sample sizes get larger
      theta(j) = bootfun(Y(g==gk(j),:));
      nk(j) = sum(g==gk(j));
      stats = bootknife(Y(g==gk(j),:),[nboot,0],bootfun,[],[],0,[],ISOCTAVE);
      SE(j) = stats.std_error;
    end
    if any(isnan(SE))
      error('evaluating bootfun on the bootknife resamples created NaN values for the standard error')
    end
    Var(j) = ((nk(j)-1)/(N-k)) * SE(j)^2;
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
    t = abs(theta(i(idx)) - theta(j(idx))) ./ sqrt(Var * (w(i(idx)) + w(j(idx))));
  else
    % Calculate Dunnett's test statistic 
    t = abs((theta - theta(ref))) ./ sqrt(Var * (w + w(ref)));
  end
  maxT = max(t);
  
end
