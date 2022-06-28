function [y, g] = unitmeans (x, clusters, nvar)

  % Helper function file required for ibootci

  % Calculate unit (cluster) means

  % Calculate number of levels of subsampling
  L = size(clusters,2);

  % Get IDs of unique clusters in lowest level
  gid = unique(clusters(:,L));
  K = numel(gid);

  % Initialize output variables
  g = zeros(K,L-1);
  y = cell(1,nvar);
  for v = 1:nvar
    y{v} = zeros(K,1);
  end

  % Calculate cluster means
  for k = 1:K

    % Find last level cluster members
    idx = find(clusters(:,L) == gid(k));

    % Compute cluster means
    for v = 1:nvar
      y{v}(k) = mean(x{v}(idx));
    end

    % Check data nesting
    if numel(unique(clusters(idx,L-1))) > 1
      error('Impossible hierarchical data structure')
    end

    % Redefine clusters
    g(k,:) = clusters(idx(1),1:L-1);

  end

end
