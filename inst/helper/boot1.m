function [T1, T2, U, idx] = boot1 (x, nboot, n, nvar, bootfun, T0, S, opt)

  % Helper function file required for ibootci

  % Extract required options structure fields
  weights = opt.weights;
  strata = opt.strata;
  blocksize = opt.blocksize;
  bandwidth = opt.bandwidth;
  R = opt.R;
  df = opt.df;
  stderr = opt.stderr;
  runmode = opt.runmode;
  paropt = opt.paropt;

  % Initialize
  B = nboot(1);
  C = nboot(2);
  N = n*B;
  T1 = zeros(1,B);
  if C>0
    T2 = zeros(C,B);
    U = zeros(1,B);
  elseif ~isempty(stderr)
    T2 = [];
    U = zeros(1,B);
  else
    T2 = [];
    U = [];
  end
  X1 = cell(1,nvar);
  if nargout < 4
    idx = zeros(n,1);
  else
    idx = zeros(n,B);
  end
  % Calculate mean and variance of the original sample
  xbar = zeros(1,nvar);
  xvar = zeros(1,nvar,1);
  for v=1:nvar
    xbar(v) = mean(x{v});
    xvar(v) = var(x{v},1);
  end

  % If applicable, prepare for stratified resampling
  if ~isempty(strata)
    % Get strata IDs
    gid = unique(strata);  % strata ID
    K = numel(gid);    % number of strata
    % Create strata matrix
    g = zeros(n,K);
    for k = 1:K
      g(:,k) = (strata == gid(k));
    end
    % Get strata sample and bootstrap sample set dimensions
    nk = sum(g).';   % strata sample sizes
    ck = cumsum(nk); % cumulative sum of strata sample sizes
    ik = [1;ck];   % strata boundaries
    Nk = nk*B;     % size of strata bootstrap sample set
    Ck = cumsum(Nk); % cumulative sum of strata bootstrap sample set sizes
  else
    ck = n;
    g = ones(n,1);
    K = 1;
    nk = n;
  end
  g = logical(g);

  % Prepare weights for resampling
  if any(diff(weights))
    if ~isempty(strata)
      % Calculate within-stratum weights
      c = zeros(n,1);
      for k = 1:K
        c = c + round(Nk(k) * g(:,k).*weights./sum(g(:,k).*weights));
        c(ik(k):ik(k+1),1) = cumsum(c(ik(k):ik(k+1),1));
        c(ik(k+1)) = Ck(k);
      end
    else
      % Calculate weights (no groups)
      c = cumsum(round(N * weights./sum(weights)));
      c(end) = N;
    end
    c = [c(1);diff(c)];
  else
    c = ones(n,1)*B;
  end

  % Perform bootstrap resampling
  try
    pool = gcp('nocreate');
  catch
    pool = [];
  end
  if paropt.UseParallel && isoctave
    % Octave parallel computing
    % Perform ordinary resampling with replacement
    if nargout > 3
      error('No bootsam when operating ibootci in parallel mode')
    end
    % Prepare resampling weights
    w = zeros(n,K);
    for k = 1:K
      w(:,k) = cumsum((g(:,k).*weights)./sum(g(:,k).*weights));
    end
    parfun = @() parboot (x, X1, nboot, n, nvar, bootfun, T0, g, S, opt, w, ck, xbar, xvar);
    bootout = cell(1,B);
    errfun = @() error('An error occurred during octave parallel implementation. Try running computations in serial.');
    bootout = parcellfun(paropt.nproc, parfun, bootout, 'UniformOutput',false,'ErrorHandler',errfun);
    T1 = cell2mat(cellfun(@(S) S.T1, bootout, 'UniformOutput', false));
    T2 = cell2mat(cellfun(@(S) S.T2.', bootout, 'UniformOutput', false));
    U = cell2mat(cellfun(@(S) S.U, bootout, 'UniformOutput', false));
  elseif (paropt.UseParallel || ~isempty(pool)) && ~isoctave
    % Matlab parallel computing
    % Perform ordinary resampling with replacement
    if nargout > 3
      error('No bootsam when operating ibootci in parallel mode')
    end
    % Prepare resampling weights
    w = zeros(n,K);
    for k = 1:K
      w(:,k) = cumsum((g(:,k).*weights)./sum(g(:,k).*weights));
    end
    parfor h = 1:B
      idx = zeros(n,1);
      X1 = cell(1,nvar);
      for i = 1:n
        k = sum(i>ck)+1;
        j = sum((rand(1) >= w(:,k)))+1;
        idx(i,1) = j;
      end
      for v = 1:nvar
        X1{v} = x{v}(idx);
      end
      % Since second bootstrap is usually much smaller, perform rapid
      % balanced resampling by a permutation algorithm
      if C>0
        [U(h), T2(:,h)] = boot2 (X1, nboot, n, nvar, bootfun, T0, g, S, opt);
      end
      if ~isempty(stderr)
        U(h) = stderr(X1{:});
      end
      if ~isempty(bandwidth)
        % Apply smoothing using a Gaussian kernel
        noise = bsxfun(@times,randn(n,nvar)*chol(R),bandwidth);
        for v = 1:nvar
          X1{v} = shrunk_smooth (X1{v}, bandwidth(v), xbar(v), xvar(v), noise(:,v));
        end
      end
      T1(h) = feval(bootfun,X1{:});
    end
  else
    % Octave or Matlab serial computing
    % Since first bootstrap is large, use a memory efficient balanced resampling algorithm
    %    Gleason, J.R. (1988) Algorithms for Balanced Bootstrap Simulations. 
    %    The American Statistician. Vol. 42, No. 4 pp. 263-266
    % If strata is provided, resampling is stratified
    % In serial mode, a random seed is specified making the Monte Carlo simulation 
    % deterministic and reproducible between Matlab and Octave
    rand('twister',rng_state);
    for h = 1:B
      for i = 1:n
        k = sum(i>ck)+1;
        j = sum((rand(1) >= cumsum((g(:,k).*c)./sum(g(:,k).*c))))+1;
        if nargout < 4
          idx(i,1) = j;
        else
          idx(i,h) = j;
        end
        c(j) = c(j)-1;
      end
      for v = 1:nvar
        if nargout < 4
          X1{v} = x{v}(idx);
        else
          X1{v} = x{v}(idx(:,h));
        end
      end
      % Since second bootstrap is usually much smaller, perform rapid
      % balanced resampling by a permutation algorithm
      if C>0
        [U(h), T2(:,h)] = boot2 (X1, nboot, n, nvar, bootfun, T0, g, S, opt);
      end
      if ~isempty(stderr)
        U(h) = stderr(X1{:});
      end
      if ~isempty(bandwidth)
        % Apply smoothing using a Gaussian kernel
        noise = bsxfun(@times,randn(n,nvar)*chol(R),bandwidth);
        for v = 1:nvar
          X1{v} = shrunk_smooth (X1{v}, bandwidth(v), xbar(v), xvar(v), noise(:,v));
        end
      end
      T1(h) = feval(bootfun,X1{:});
    end
  end

end
