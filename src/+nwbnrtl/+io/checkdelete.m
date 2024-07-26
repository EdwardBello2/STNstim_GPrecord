function [TF] = checkdelete(fullpath)
% Checks if a given NWB file exists at the full path specified by user. If
% so, deletes it. Outputs TF. True: file found and deleted. False: file not
% found. Only works for detecting/deleting NWB files. 

TF = false;

% append nwb extension if not specified
[filepath, name, ext] = fileparts(fullpath);
if isempty(ext)
    fullpath = [filepath name '.nwb'];
    
end

% delete the file if it exists
if isfile(fullpath)
    delete(fullpath);
    TF = true;
    
end


end