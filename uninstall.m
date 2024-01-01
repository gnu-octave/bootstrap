% Basic uninstall script for local installation
% 

% Get inst directory
inst_dir = fullfile (fileparts (mfilename ('fullpath')), 'inst');

% Check if running in Octave (else assume Matlab)
info = ver; 
isoctave = any (ismember ({info.Name}, 'Octave'));

if isoctave
  % Uninstall for Octave
  rmpath (inst_dir)
  octaverc = '~/.octaverc';
  if (exist (octaverc,'file'))
    [fid, msg] = fopen (octaverc, 'r+t');
    S = (fread (fid, '*char')).';
    fclose(fid);
    [fid, msg] = fopen (octaverc, 'wt');
  else
    error('~/.octaverc does not exist');
  end
  comment = regexptranslate ('escape', '% Load statistics-resampling package');
  S = regexprep(S,['\r\n\r\n',comment],'');
  S = regexprep(S,strcat('\r\n',...
                  regexptranslate ('escape', strcat('addpath (''', inst_dir, ''', ''-end'')'))),'');
  fseek (fid, 0);
  fputs (fid, S);
  fclose (fid);
else
  % Assumming uninstall for Matlab instead
  rmpath (inst_dir)
  if exist('savepath')
    savepath
  else
    % backwards compatibility
    path2rc;
  end
end


% Notify user that uninstall is complete
disp ('This statistics-resampling package has been uninstalled from this location')

% Clean up
clear info isoctave S comment octaverc fid msg inst_dir
