% Check if running in Octave (else assume Matlab)
info = ver; 
isoctave = any (ismember ({info.Name}, 'Octave'));

% Create cell array of paths to add
arch = computer('arch');
[comp, maxsize, endian] = computer ();
binary = true;
switch endian
  case 'L' 
    if isoctave
      arch_idx = regexpi (arch, {'mingw32-i686', 'mingw32-x86_64', '^darwin.*x86_64$', 'gnu-linux-x86_64'});
      binary_paths = {'..\bin\octave\win\win32\','..\bin\octave\win\win64\','../bin/octave/mac/maci64/', '../bin/octave/linux/glnxa64/'};
    else
      arch_idx = regexpi (arch, {'win32', 'win64', 'maci64', 'glnxa64'});
      binary_paths = {'..\bin\matlab\win\win32\','..\bin\matlab\win\win64\','../bin/matlab/mac/maci64/', '../bin/matlab/linux/glnxa64/'};
    end
    if ~all(cellfun(@isempty,arch_idx))
      copyfile (sprintf ('%s%s%s', binary_paths{~cellfun(@isempty,arch_idx)}, 'boot', mexext), sprintf ('%s%s%s%s', 'inst', filesep, 'boot', mexext), 'f');
      copyfile (sprintf ('%s%s%s', binary_paths{~cellfun(@isempty,arch_idx)}, 'boot', mexext), sprintf ('%s%s%s%s%s%s', 'inst', filesep, 'param', filesep, 'boot', mexext), 'f');
    else
      binary = false;
    end
  case 'B'
    binary = false;
end
dirlist = {'helper', 'param', 'legacy', sprintf ('legacy%shelper', filesep)};
dirname = fileparts (mfilename ('fullpath'));

% Attemt to compile binaries from source code automatically if no suitable binaries can be found
if ~binary
  disp('No precombined binaries available for this architecture. Attempting to compile from source code...');
  if isoctave
    try
      mkoctfile --mex --output ./inst/boot ./src/boot.cpp
    catch
      err = lasterror();
      disp(err.message);
      warning ('Could not compile boot.%s. Falling back to the (slower) boot.m file.',mexext)
    end
    path_to_smoothmedian = sprintf ('./inst/param/smoothmedian.%s',mexext);
    try
      mkoctfile --mex --output ./inst/param/smoothmedian ./src/smoothmedian.cpp
    catch
      err = lasterror();
      disp(err.message);
      warning ('Could not compile smoothmedian.%s. Falling back to the (slower) smoothmedian.m file.',mexext)
    end
  else
    try  
      mex -setup c++
    catch
      err = lasterror();
      disp(err.message);
    end
    try
      mex -compatibleArrayDims -output ./bin/boot ./src/boot.cpp
    catch
      err = lasterror();
      disp(err.message);
      warning ('Could not compile boot.%s. Falling back to the (slower) boot.m file.',mexext)
    end
    try
      mex -compatibleArrayDims -output ./bin/smoothmedian ./src/smoothmedian.cpp
    catch
      err = lasterror();
      disp(err.message);
      warning ('Could not compile smoothmedian.%s. Falling back to the (slower) smoothmedian.m file.',mexext)
    end
  end
end