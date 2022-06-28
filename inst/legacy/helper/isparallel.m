function retval = isparallel

  % Helper function file required for ibootci

  % Check we have parallel computing capabilities
  str = '^parallel';
  if isoctave
    software = pkg('list');
    names = cellfun(@(S) S.name, software, 'UniformOutput', false);
    status = cellfun(@(S) S.loaded, software, 'UniformOutput', false);
    index = find(~cellfun(@isempty,regexpi(names,str)));
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
    software = {software.Name};
    if isempty(cell2mat(regexpi(software,str)))
      retval = false;
    else
      retval = true;
    end
  end
