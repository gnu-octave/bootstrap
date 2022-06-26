function bootout = parboot (idx, x, X1, nboot, n, nvar, bootfun, T0, g, S, opt, xbar, xvar)

  % Helper function file required for Octave Parallel Computing
  % capabilities of ibootci.

  % Extract required options structure fields
  blocksize = opt.blocksize;
  bandwidth = opt.bandwidth;
  R = opt.R;
  stderr = opt.stderr;

  % Resampling (with weights if applicable)
  for v = 1:nvar
    X1{v} = x{v}(idx);
  end

  % Since second bootstrap is usually much smaller, perform rapid
  % balanced resampling by a permutation algorithm
  C = nboot(2);
  if C>0
    [U, T2] = boot2 (X1, nboot, n, nvar, bootfun, T0, g, S, opt);
  else
    T2 = [];
    U = [];
  end
  if ~isempty(stderr)
    U = stderr(X1{:});
  end
  if ~isempty(bandwidth)
    % Apply smoothing using a Gaussian kernel
    noise = bsxfun(@times,randn(n,nvar)*chol(R),bandwidth);
    for v = 1:nvar
      X1{v} = shrunk_smooth (X1{v}, bandwidth(v), xbar(v), xvar(v), noise(:,v));
    end
  end
  T1 = feval(bootfun,X1{:});

  % Prepare output structure
  bootout.T1 = T1;
  bootout.T2 = T2;
  bootout.U = U;

end
