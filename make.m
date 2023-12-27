% Script to make select the appropriate precompiled binaries for installation or make the mex files from source 

% Check if running in Octave (else assume Matlab)
info = ver; 
isoctave = any (ismember ({info.Name}, 'Octave'));

% Initialise error flag
errflag = false;

% Create cell array of paths to add
arch = computer('arch');
[comp, maxsize, endian] = computer ();
try
  switch endian
    case 'L' 
      if isoctave
        arch_idx = regexpi (arch, {'mingw32-i686',...
                                   'mingw32-x86_64',... 
                                   '^darwin.*x86_64$',... 
                                   'gnu-linux-x86_64'});
        binary_paths = {'precompiled_binaries\octave\win\win32\',...
                        'precompiled_binaries\octave\win\win64\',...
                        'precompiled_binaries/octave/mac/maci64/',... 
                        'precompiled_binaries/octave/linux/glnxa64/'};
      else
        arch_idx = regexpi (arch, {'win32',... 
                                   'win64',... 
                                   'maci64',... 
                                   'glnxa64'});
        binary_paths = {'precompiled_binaries\matlab\win\win32\',...
                        'precompiled_binaries\matlab\win\win64\',...
                        'precompiled_binaries/matlab/mac/maci64/',... 
                        'precompiled_binaries/matlab/linux/glnxa64/'};
      end
      if ( ~ all (cellfun (@isempty,arch_idx)))
        retval = '0';
        while ~ismember (retval,{'1','2',''})
          retval = input (sprintf(['Potentially compatible precompiled mex files found. ', ...
                                   'Please select an option:\n', ...
                                   ' 1) Use the precompiled mex files.\n', ...
                                   ' 2) Compile new mex files from source (default).\n', ...
                                   'Answer (1 or 2): ']), 's');
        end
        if isempty(retval)
          retval = '2';
        end
        if (strcmpi (retval,'1'))
          copyfile (sprintf (repmat ('%s',1,3), binary_paths{~cellfun (@isempty, arch_idx)}, 'boot.', mexext),...
                    sprintf (repmat ('%s',1,6), '.', filesep, 'inst', filesep, 'boot.', mexext), 'f');
          copyfile (sprintf (repmat ('%s',1,3), binary_paths{~cellfun (@isempty, arch_idx)}, 'smoothmedian.', mexext),... 
                    sprintf (repmat ('%s',1,8), '.', filesep, 'inst', filesep, 'param', filesep, 'smoothmedian.', mexext), 'f');
          % Perform basic tests
          % If tests fail, try compiling source code instead
          clear boot smoothmedian
          cd inst
          boot (1, 1);
          cd param
          smoothmedian (1);
          binary = true;
          cd ../..
        else
          error ('Break from try-catch statement')
        end
      else
        binary = false;
      end
  case 'B'
      binary = false;
  end
catch
  binary = false;
end
dirlist = {'helper', 'param', 'legacy', sprintf('legacy%shelper', filesep)};
dirname = fileparts (mfilename ('fullpath'));

% Attemt to compile binaries from source code automatically if no suitable binaries can be found
if binary
  fprintf('The following suitable binaries were detected and copied over to the inst directory: ');
  fprintf('\n%s%s%s%s%s', '.', filesep, binary_paths{~cellfun(@isempty,arch_idx)}, 'boot.', mexext);
  fprintf('\n%s%s%s%s%s', '.', filesep, binary_paths{~cellfun(@isempty,arch_idx)}, 'smoothmedian.', mexext);
else
  disp('Either you chose to compile from source, or no binaries are suitable.');
  disp('Attempting to compile the source code...');
  if isoctave
    try
      mkoctfile --mex --output ./inst/boot ./src/boot.cpp
    catch
      errflag = true;
      err = lasterror();
      disp(err.message);
      warning ('Could not compile boot.%s. Falling back to the (slower) boot.m file.',mexext)
    end
    path_to_smoothmedian = sprintf('./inst/param/smoothmedian.%s',mexext);
    try
      mkoctfile --mex --output ./inst/param/smoothmedian ./src/smoothmedian.cpp
    catch
      errflag = true;
      err = lasterror();
      disp(err.message);
      warning ('Could not compile smoothmedian.%s. Falling back to the (slower) smoothmedian.m file.',mexext)
    end
    try
      mkoctfile --mex --output ./inst/boot ./src/boot.cpp
    catch
      errflag = true;
      err = lasterror();
      disp(err.message);
      warning ('Could not compile boot.%s. Falling back to the (slower) boot.m file.',mexext)
    end
  else
    try  
      mex -setup c++
    catch
      errflag = true;
      err = lasterror();
      disp(err.message);
    end
    try
      mex -compatibleArrayDims -output ./inst/boot ./src/boot.cpp
    catch
      errflag = true;
      err = lasterror();
      disp(err.message);
      warning ('Could not compile boot.%s. Falling back to the (slower) boot.m file.',mexext)
    end
    try
      mex -compatibleArrayDims -output ./inst/param/smoothmedian ./src/smoothmedian.cpp
    catch
      errflag = true;
      err = lasterror();
      disp(err.message);
      warning ('Could not compile smoothmedian.%s. Falling back to the (slower) smoothmedian.m file.',mexext)
    end
  end
end
if errflag
  fprintf('\nmake completed with errors. Please review the details in the errors in the above output. \n')
else
  fprintf('\n''make'' completed. Please now run the ''install'' command. \n')
end


clear arch arch_idx binary binary_paths comp dirlist dirname endian info isoctave maxsize
