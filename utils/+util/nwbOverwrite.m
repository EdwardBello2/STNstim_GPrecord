function nwbUpd = nwbOverwrite(nwb, filename)
% Overwrite original nwb file with current updated nwb object.
% Intended for effectively adding info to current nwb file. This function
% is specifically for "updating" a file without changing its name; if you 
% want to create a new unique nwb file, just use nwbExport.

% Details: First creates an updated "copy" of the original file with any 
% alterations, saving it with the same name + "_update". Then deletes the
% original file and renames the new one without the "_update" suffix. In
% effect, it *looks* like the file just got changes saved to it. In reality
% a brand new file was created and then renamed. 

% Assert that filename already exists
if ~isfile(filename)
    error('Specified file does not exist. Remember to include extention')
    
end

% Create "new" nwb file
[filepath, name, ext] = fileparts(filename);
filenameUpd = [filepath '\' name '_update' ext];
nwbExport(nwb, filenameUpd)


% delete "original" nwb file
delete(filename)


% rename "new" nwb file 
movefile(filenameUpd, filename)

nwbUpd = nwbRead(filename);

end