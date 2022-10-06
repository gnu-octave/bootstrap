function stats = bootblock(data,blocksize,nboot,alpha)

    sz = size(data);
    n = sz(1)
 
    % Prepare for block resampling (if applicable)
    if strcmpi(blocksize,'auto')
      % Set block length
      blocksize = max(round(n^(1/3)),2);   % in the order of n^(1/3)
    end
    data = split_blocks(data,n,blocksize)


end

function y = split_blocks (x, l1, l2)

  % Helper function file required for ibootci

  % Calculate data and block dimensions
  n = size(x,1);
  nvar = numel(x);

  % Create a matrix of circular, overlapping blocks
  % Ref: Politis and Romano (1991) Technical report No. 370
  y = zeros(n,l2);
  l = l1;
  i = 1;
  while l > 0
    temp = cat(1, x(((i-1)*l1+1:min(n,i*l1)),1), ...
                  x((i-1)*l1+1:min(n,(i-1)*l1+l2),1));
    for k = 1:l
      y((i-1)*l1+k,:) = temp(k:k+l2-1);
    end
    l = min(l1,n-(i*l1));
    i = i + 1;
  end
  y

end

function y = cat_blocks (nvar, varargin)

  % Helper function required for block bootstrap in ibootci 

  % Get data dimensions
  x = (varargin);
  N = numel(x);
  l = N/nvar;
  [n, reps] = size(x{1});

  % Concatenate blocks
  y = cell(1,nvar);
  y = zeros(N,reps);
  for i = 1:l
    y(i:l:n*l,:) = x{(v-1)*l+i};
  end
  y = y(1:n,:);

end
