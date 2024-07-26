function [fullpath] = remextension(fullpath)
% function to remove the extension part of a string describing a file
[pn, name, ~] = fileparts(fullpath);

if isempty(pn)
    fullpath = name;
    
else
    fullpath = [pn '\' name];

end

end