% Basic uninstall script 
% 

copyfile ('PKG_ADD','PKG_ADD.m'); % copy to m-file for backwards compatibility
copyfile ('PKG_DEL','PKG_DEL.m'); % copy to m-file for backwards compatibility
run (fullfile(pwd,'PKG_ADD.m'));

% Check if running in Octave (else assume Matlab)
info = ver; 
isoctave = any (ismember ({info.Name}, 'Octave'));

if isoctave
  % Uninstall for Octave
  dirlist = cell(5,1); % dir list needs to be in decreasing order of length
  dirlist{1} = fullfile (pwd,'inst','legacy','helper');
  dirlist{2} = fullfile (pwd,'inst','legacy');
  dirlist{3} = fullfile (pwd,'inst','helper');
  dirlist{4} = fullfile (pwd,'inst','param');
  dirlist{5} = fullfile (pwd,'inst','');
  n = numel (dirlist);
  octaverc = '~/.octaverc';
  if exist(octaverc,'file')
    [fid, msg] = fopen (octaverc, 'r+t');
    S = (fread (fid, '*char')).';
    fclose(fid);
    [fid, msg] = fopen (octaverc, 'wt');
  else
    error('~/.octaverc does not exist');
  end
  comment = regexptranslate ('escape', '% Load statistics-bootstrap package');
  S = regexprep(S,['\r\n\r\n',comment],'');
  for i=1:n
    S = regexprep(S,strcat('\r\n',...
                  regexptranslate ('escape', strcat('addpath (''',dirlist{i},''', ''-end'')'))),'');
  end
  fseek (fid, 0);
  fputs (fid, S);
  fclose (fid);
  % Unload package from current session
  run (fullfile(pwd,'PKG_DEL.m'));
else
  % Assumming uninstall for Matlab instead
  run (fullfile(pwd,'PKG_DEL'));
  if exist('savepath')
    savepath
  else
    % backwards compatibility
    path2rc;
  end
end

% Delete MEX files
path_to_boot = sprintf ('./inst/boot.%s',mexext);
if exist (path_to_boot)
  delete (path_to_boot);
  clear boot
end
path_to_smoothmedian = sprintf ('./inst/param/smoothmedian.%s',mexext);
if exist (path_to_smoothmedian)
  delete (path_to_smoothmedian);
  clear smoothmedian
end


% Notify user that uninstall is complete
disp ('This statistics-bootstrap package has been uninstalled from this location')

% Clean up
clear dirlist S comment i ii octaverc fid n msg
delete ('PKG_ADD.m');
delete ('PKG_DEL.m');
