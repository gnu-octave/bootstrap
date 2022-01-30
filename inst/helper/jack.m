function [SE, T, U] = jack (x, func, paropt, opt)

  % Helper function file required for ibootci

  % Jackknife

  % Check what type of variable x is
  if ~iscell(x)
    x = num2cell(x,1); % convert to cell array
    matflag = true;
    opt.matflag = matflag;
  else 
    matflag = false;
    opt.matflag = matflag;
  end

  % Get basic info about data
  nvar = size(x,2);
  [m, n] = size(x{1});

  if nargin < 2
    error('Invalid number of input arguments');
  end
  if nargin < 3 || isempty(paropt)
    paropt = struct;
    paropt.UseParallel = false;
    paropt.nproc = 1;
  end
  if nargin < 4 || isempty(opt)
    opt = struct;
    opt.weights = [];
    opt.blocksize = [];
    opt.clusters = [];
  end
  if nargout > 3
    error('Invalid number of output arguments');
  end

  % Prepare data for resampling
  if nargin > 3
    if isfield(opt,'blocksize') && ~isempty(opt.blocksize)
      % Prepare for Block-Jackknife
      % Pad data (circular)
      % Overlapping blocks of size l will be primary sampling unit
      l = opt.blocksize;
      for v = 1:nvar
        x{v} = num2cell(cat(1,x{v},x{v}(1:l-1)));  % for circular blocks
      end
    elseif isfield(opt,'clusters') && ~isempty(opt.clusters)
      % Prepare Cluster-Jackknife
      % Set whole clusters as primary sampling units
      [~, idx] = sort(opt.clusters);
      [~, ~, m, g] = sse_calc (x, opt.clusters, nvar);
      l = 1;
      for v = 1:nvar
        x{v} = x{v}(idx);
        x{v} = mat2cell(x{v},sum(g));
      end
    else
      l = 1;
      % Prepare for Jackknife
      % Each original data value is a primary sampling unit
      for v = 1:nvar
        x{v} = num2cell(x{v});
      end
    end
    if isfield(opt,'weights') && ~isempty(opt.weights)
      % Prepare weights
      if any(diff(opt.weights)) && isempty(opt.clusters)
        w = opt.weights/sum(opt.weights);
      else
        w = 1/m .* ones(m,1);
      end
    else 
      w = 1/m .* ones(m,1);
    end
  else
    l = 1;
    % Prepare for Jackknife
    % Each original data value is a primary sampling unit
    for v = 1:nvar
      x{v} = num2cell(x{v});
    end
    w = 1/m .* ones(m,1);
  end

  % Perform 'leave one out' procedure and calculate the variance(s)
  % of the test statistic.
  T = zeros(m,n);
  try
    pool = gcp('nocreate');
  catch
    pool = [];
  end
  if paropt.UseParallel && isoctave

    % Octave parallel computing
    i = [1:m].';
    parfun = @(i) parjack (x, func, nvar, m, i, l);
    T = pararrayfun(paropt.nproc, parfun, i, 'UniformOutput',true);

  elseif (paropt.UseParallel || ~isempty(pool)) && ~isoctave

    % Matlab parallel computing
    parfor i = 1:m
      M = cell(1,nvar);
      for v = 1:nvar
        M{v} = x{v};
        M{v}(i:(i+l-1),:) = [];
        M{v} = cell2mat(M{v});
      end
      if opt.matflag > 0
        T(i,:) = feval(func, cell2mat(M));
      else
        T(i,:) = feval(func, M{:});
      end
      M = [];
    end

  else

    % Octave or Matlab serial computing
    for i = 1:m
      M = cell(1,nvar);
      for v = 1:nvar
        M{v} = x{v};
        M{v}(i:(i+l-1),:) = [];
        M{v} = cell2mat(M{v});
      end
      if opt.matflag > 0
        T(i,:) = feval(func, cell2mat(M));
      else
        T(i,:) = feval(func, M{:});
      end
      M = [];
    end

  end
  
  % Calculate influence function and sampling variance for func
  % from the Jackknife leave-one-out statistics
  Tori = sum(w .* T, 1);
  Tori = Tori(ones(m, 1), :);
  Ti = (m * Tori - (m - l) * T) / l;
  % Note that (Ti - Tori) / (m - l) == (Tori - T)
  U = (m - l) * (Tori - T);
  Var = sum(w .* (U.^2)/(m - 1), 1);

  % Calculate standard error(s) of the functional parameter
  SE = sqrt(Var);

end
