% Basic script for local installation.
% Make sure you run the 'make' function first for optimum performance.
%

% Add inst directory to path
inst_dir = fullfile (fileparts (mfilename ('fullpath')), 'inst');
addpath (inst_dir)

% Check if running in Octave (else assume Matlab)
info = ver; 
isoctave = any (ismember ({info.Name}, 'Octave'));

% Save the newly added paths so that they will be loaded each time we start Octave or Matlab
if isoctave
  % Install for Octave
  octaverc = '~/.octaverc';
  if exist(octaverc,'file')
    [fid, msg] = fopen (octaverc, 'r+t');
  else
    [fid, msg] = fopen (octaverc, 'w+t');
  end 
  S = (fread (fid, '*char')).';
  comment = sprintf ('\r\n\r\n%s', '% Load statistics-resampling package');
  S = strcat (S, comment);
  S = strcat (S, sprintf ('\r\naddpath (''%s'', ''-end'')', inst_dir));
  fseek (fid, 0);
  fputs (fid, S);
  fclose (fid);
else
  % Assuming install for Matlab instead
  if exist('savepath')
    savepath;
  else
    % backwards compatibility
    path2rc;
  end
end

% Notify user that installation is complete
d = struct;
d.dir = inst_dir;
post_install (d);

% Check if the user has run the make script
inst_files = dir (inst_dir);
if (all (arrayfun (@(name) ismember (sprintf ('%s.%s', name{:}, mexext), ...
                                 {inst_files.name}), {'boot', 'smoothmedian'})))
  try
    boot (1, 1);
    smoothmedian (1);
    make_done = true;
  catch
    make_done = false;
  end
else
  make_done = false;
end
if (~ make_done)
  warning ('For optimal performance, run the ''make'' command to copy or compile the appropriate binaries')
end

% Clean up
clear info isoctave S comment octaverc fid msg inst_dir inst_files d make_done

