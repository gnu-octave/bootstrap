% Basic uninstall script 
% 

copyfile ('PKG_ADD','PKG_ADD.m'); % copy to m-file for backwards compatibility
copyfile ('PKG_DEL','PKG_DEL.m'); % copy to m-file for backwards compatibility
run (fullfile(pwd,'PKG_ADD.m'));
if isoctave
  % Uninstall for Octave
  dirlist = cell(3,1); % dir list needs to be in decreasing order of length
  dirlist{1} = fullfile (pwd,'inst','helper');
  dirlist{2} = fullfile (pwd,'inst','param');
  dirlist{3} = fullfile (pwd,'inst','');
  n = numel(dirlist);
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
  if exist ('./inst/boot.oct')
    delete ('./inst/boot.oct');
  end
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

% Unload package from current session
run (fullfile(pwd,'PKG_DEL.m'));

% Notify user that uninstall is complete
disp ('This statistics-bootstrap package has been uninstalled from this location')

% Clean up
clear dirlist S comment i ii octaverc fid n msg
delete ('PKG_ADD.m');
delete ('PKG_DEL.m');
