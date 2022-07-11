% Basic script for local installation
% 

% Run make function to copy over the appropriate binaries to the inst directory, or compile the binaries from the source code
make;

% Add directories to path
copyfile ('PKG_ADD','PKG_ADD.m');
run (fullfile(pwd,'PKG_ADD.m'));
dirlist = cell(5,1); % dir list needs to be in increasing order of length
dirlist{1} = fullfile (pwd,'inst','');
dirlist{2} = fullfile (pwd,'inst','param');
dirlist{3} = fullfile (pwd,'inst','helper');
dirlist{4} = fullfile (pwd,'inst','legacy');
dirlist{5} = fullfile (pwd,'inst','legacy','helper');
n = numel (dirlist);

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
  comment = sprintf ('\r\n\r\n%s', '% Load statistics-bootstrap package');
  if isempty(strfind(S,comment))
    S = strcat (S, comment);
  end
  for i=1:n
    if isempty(strfind(S,dirlist{i}))
      S = strcat (S, sprintf ('\r\naddpath (''%s'', ''-end'')', dirlist{i}));
    end
  end
  fseek (fid, 0);
  fputs (fid, S);
  fclose (fid);
else
  % Assumming install for Matlab instead
  if exist('savepath')
    savepath;
  else
    % backwards compatibility
    path2rc;
  end
end

% Notify user that installation is complete
disp ('The statistics-bootstrap package has been installed at the current location ')

% Clean up
clear info isoctave dirlist S comment i ii octaverc fid n msg
delete ('PKG_ADD.m');
