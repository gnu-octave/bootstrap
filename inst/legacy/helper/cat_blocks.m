function y = cat_blocks (nvar, varargin)

  % Helper function required for block bootstrap in ibootci 

  % Get data dimensions
  x = (varargin);
  N = numel(x);
  l = N/nvar;
  [n, reps] = size(x{1});

  % Concatenate blocks
  y = cell(1,nvar);
  for v = 1:nvar
    y{v} = zeros(N,reps);
    for i = 1:l
      y{v}(i:l:n*l,:) = x{(v-1)*l+i};
    end
    y{v} = y{v}(1:n,:);
  end

end
