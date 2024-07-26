function cachedDataPath = assertCacheDir(cfg, scriptName)
% Enter the name of the currently-executing m-file; if a folder for it does
% not yet exist in the cachedData directory, create it. 

% global PROJROOTPATH

% [scriptDirectoryFullPath, scriptName] = fileparts(mfilename('fullpath'));

% First make sure that this script has a folder within the project folder's
% intermediate data section
% rootDir = assertEndSlash(cfg.projRootDirFullpath);
cacheDir = assertEndSlash(cfg.cacheDirFullpath);

cachedDataPath = [cacheDir scriptName];
if ~exist(cachedDataPath, 'dir')
    mkdir(cachedDataPath)
    
end

cachedDataPath = [cachedDataPath '\'];

end
