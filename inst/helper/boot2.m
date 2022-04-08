function [U, T2] = boot2 (X1, nboot, n, nvar, bootfun, T0, g, S, opt)

  % Helper function file required for ibootci

  % Extract required options structure fields
  blocksize = opt.blocksize;
  runmode = opt.runmode;

  % Note that weights are not implemented here with iterated bootstrap

  % Prepare for block resampling (if applicable)
  if ~isempty(blocksize)
    x1 = cat_blocks(S.nvar,X1{:});
    blocksize = max(fix(blocksize/2),1);
    X1 = split_blocks(x1,opt.blocksize,blocksize);
    nvar = S.nvar * blocksize;
    g = ones(n,1);
  end

  % Initialize
  C = nboot(2);

  % If applicable, prepare for stratified resampling
  K = size(g,2);    % number of strata
  nk = sum(g).';    % strata sample sizes
  ck = cumsum(nk);  % cumulative sum of strata sample sizes
  ik = [1;ck+1];    % strata boundaries (different definition to boot1)
  Nk = nk*C;        % size of strata bootstrap sample set

  % Rapid balanced resampling by permutation
  % If strata is provided, resampling is stratified
  idx = zeros(n,C);
  for k = 1:K
    tmp = (1:nk(k))'*ones(1,C);
    [~,rperm] = sort(rand(1,Nk(k)));
    tmp = tmp(reshape(rperm,nk(k),C));
    tmp = tmp + ik(k) - 1;
    idx(ik(k): ik(k+1)-1,:) = tmp;
  end
  X2 = cell(1,nvar);
  for v = 1:nvar
    X2{v} = X1{v}(idx);
  end
  switch lower(runmode)
    case {'fast'}
      % Vectorized calculation of second bootstrap statistics
      T2 = feval(bootfun,X2{:});
    case {'slow'}
      % Calculation of second bootstrap statistics using a loop
      T2 = zeros(1,C);
      for i=1:C
        x2 = cellfun(@(X2)X2(:,i),X2,'UniformOutput',false);
        T2(i) = feval(bootfun,x2{:});
      end
  end
  U = interp_boot2(T2,T0,C);

end
