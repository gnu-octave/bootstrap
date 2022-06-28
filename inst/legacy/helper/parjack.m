function T = parjack (x, func, nvar, m, i, l)

  % Helper function file required for Octave Parallel Computing
  % capabilities of jack.

  M = cell(1,nvar);
  for v = 1:nvar
    M{v} = x{v};
    M{v}(i:(i+l-1),:) = [];
    M{v} = cell2mat(M{v});
  end
  T = feval(func, M{:});
  M = [];

end
