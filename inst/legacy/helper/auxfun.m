function T = auxfun (bootfun, nvar, varargin)

  % Helper function required for block bootstrap in ibootci
  X = varargin{1};
  Y = cat_blocks(nvar,X{:});
  T = bootfun(Y{:});

end
