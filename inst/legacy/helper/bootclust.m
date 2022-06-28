function T = bootclust (bootfun, K, g, runmode, mu, varargin)

  % Helper function required for cluster bootstrap in ibootci

  % Two-stage nonparametric bootstrap sampling with shrinkage
  % correction for clustered data [1-4].
  %
  % By resampling residuals, this bootstrap method can be used when
  % cluster sizes are unequal. However, cluster samples are assumed
  % to be taken from populations with equal variance. Not compatible
  % with bootstrap-t or bootstrap iteration.
  %
  % References:
  %  [1] Davison and Hinkley (1997) Bootstrap Methods and their
  %       application. Chapter 3: pg 97-100
  %  [2] Ng, Grieve and Carpenter (2013) The Stata Journal.
  %       13(1): 141-164
  %  [3] Gomes et al. (2012) Medical Decision Making. 32(2): 350-361
  %  [4] Gomes et al. (2012) Health Econ. 21(9):1101-18

  % Calculate data dimensions
  Z = varargin{1};
  nvar = numel(Z);
  [n,reps] = size(Z{1});

  % Preallocate arrays
  bootmu = cell(1,nvar);
  X = cell(1,nvar);
  for v = 1:nvar
    X{v} = zeros(n,reps);
  end

  % Ordinary resampling with replacement of cluster means
  idx = ceil(K.*rand(K,reps));   % For compatibility with R2007
  for v = 1:nvar
    bootmu{v} = mu{v}(idx);
  end

  % Combine residuals with resampled cluster means
  for v = 1:nvar
    for k = 1:K
      X{v}(g(:,k),:) = bsxfun(@plus, Z{v}(g(:,k),:), bootmu{v}(k,:));
    end
  end

  % Calculate bootstrap statistic(s)
  switch lower(runmode)
    case {'fast'}
      T = feval(bootfun,X{:});
    case {'slow'}
      T = zeros(1,reps);
      for i = 1:reps
        x = cellfun(@(X) X(:,i),X,'UniformOutput',false);
        T(i) = feval(bootfun,x{:});
      end
  end

end
