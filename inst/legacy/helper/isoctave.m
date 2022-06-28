function retval = isoctave

  % Helper function file required for ibootci

  % Check we are running within Octave
  software = ver;
  software = {software.Name};
  if any(strcmpi(software,'OCTAVE'))
    retval = true;
  else
    retval = false;
  end

end
