function retval = isparallel

  % Helper function file required for ibootci

  % Check we have parallel computing capabilities
  if isoctave
    pat = '^parallel';
    software = pkg('list');
    names = cellfun(@(S) S.name, software, 'UniformOutput', false);
    status = cellfun(@(S) S.loaded, software, 'UniformOutput', false);
    index = find(~cellfun(@isempty,regexpi(names,pat)));
    if ~isempty(index)
      if logical(status{index})
        retval = true;
      else
        retval = false;
      end
    else
      retval = false;
    end
  else
    software = ver;
    if ismember ('Parallel Computing Toolbox', {software.Name})
      retval = true;
    else
      retval = false;
    end
  end
