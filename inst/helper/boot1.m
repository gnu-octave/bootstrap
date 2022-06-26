function [T1, T2, U, idx] = boot1 (x, nboot, n, nvar, bootfun, T0, S, opt)

  % Helper function file required for ibootci

  % Extract required options structure fields
  weights = opt.weights;
  strata = opt.strata;
  blocksize = opt.blocksize;
  bandwidth = opt.bandwidth;
  R = opt.R;
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
    K = numel(gid);        % number of strata
    % Create strata matrix
    g = zeros(n,K);
    for k = 1:K
      g(:,k) = (strata == gid(k));
    end
    % Get strata sample and bootstrap sample set dimensions
    nk = sum(g).';         % strata sample sizes
    ck = cumsum(nk);       % cumulative sum of strata sample sizes
    ik = [1;ck];           % strata boundaries
    Nk = nk*B;             % size of strata bootstrap sample set
    Ck = cumsum(Nk);       % cumulative sum of strata bootstrap sample set sizes
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

  % Perform balanced bootknife resampling
  if ~isempty (strata)
    idx = zeros (n, B, 'int16');
    for k = 1:K
      idx(g(:, k),:) = boot (nk(k), B, false, c(g(:,k)));
      rows = find (g(:, k));
      idx(g(:, k),:) = rows(idx(g(:, k), :));
    end
  else
    idx = boot (n, B, false, c);
  end
  try
    pool = gcp('nocreate');
  catch
    pool = [];
  end
  if paropt.UseParallel && isoctave
    % Octave parallel computing
    parfun = @(idx) parboot (idx, x, X1, nboot, n, nvar, bootfun, T0, g, S, opt, xbar, xvar);
    errfun = @() error('An error occurred during octave parallel implementation. Try running computations in serial.');
    bootout = parcellfun(paropt.nproc, parfun, num2cell(idx,1), 'UniformOutput', false);
    T1 = cell2mat(cellfun(@(S) S.T1, bootout, 'UniformOutput', false));
    T2 = cell2mat(cellfun(@(S) S.T2.', bootout, 'UniformOutput', false));
    U = cell2mat(cellfun(@(S) S.U, bootout, 'UniformOutput', false));
  elseif (paropt.UseParallel || ~isempty(pool)) && ~isoctave
    % Matlab parallel computing
    parfor h = 1:B
      X1 = cell(1,nvar);
      for v = 1:nvar
        X1{v} = x{v}(idx(:,h));
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
    for h = 1:B
      X1{v} = x{v}(idx(:,h));
      % If applicable, perform second bootstrap 
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
