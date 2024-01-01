% Basic script for local installation
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
disp ('The statistics-resampling package has been installed at the current location ')

% Clean up
clear info isoctave S comment octaverc fid msg inst_dir

