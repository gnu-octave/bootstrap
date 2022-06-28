function y = split_blocks (x, l1, l2)

  % Helper function file required for ibootci

  % Calculate data and block dimensions
  n = size(x{1},1);
  nvar = numel(x);

  % Create a matrix of circular, overlapping blocks
  % Ref: Politis and Romano (1991) Technical report No. 370
  y = cell(1,nvar);
  for v = 1:nvar
    y{v} = zeros(n,l2);
    l = l1;
    i = 1;
    while l > 0
      temp = cat(1, x{v}(((i-1)*l1+1:min(n,i*l1)),1), ...
                    x{v}((i-1)*l1+1:min(n,(i-1)*l1+l2),1));
      for k = 1:l
        y{v}((i-1)*l1+k,:) = temp(k:k+l2-1);
      end
      l = min(l1,n-(i*l1));
      i = i + 1;
    end
  end
  y = cell2mat(y);
  y = num2cell(y,1);

end
