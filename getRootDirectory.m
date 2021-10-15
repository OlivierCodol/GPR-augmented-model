function dir = getRootDirectory

script = evalin('caller','matlab.desktop.editor.getActiveFilename');
idx = strfind(script,'\');
dir = script(1:(idx(end)-1));

end

