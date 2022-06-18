% Basic script for local installation
% 

copyfile ('PKG_ADD','PKG_ADD.m');
run (fullfile(pwd,'PKG_ADD.m'));
dirlist = cell(3,1); % dir list needs to be in increasing order of length
dirlist{1} = fullfile (pwd,'inst','');
dirlist{2} = fullfile (pwd,'inst','param');
dirlist{3} = fullfile (pwd,'inst','helper');
n = numel (dirlist);
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
  try
    mkoctfile --output ./inst/boot.oct ./src/boot.cc
  catch
    warning ('Could not compile boot.oct. Falling back to the (slower) boot.m file.')
  end
  try
    mkoctfile --output ./inst/smoothmedian.oct ./src/smoothmedian.cc
  catch
    warning ('Could not compile smoothmedian.oct. Falling back to the (slower) smoothmedian.m file.')
  end
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
clear dirlist S comment i ii octaverc fid n msg
delete ('PKG_ADD.m');
