% script to create a SARGE-compatible .mat file for DBS artifact detection

% script to change .mat "tdt" format of data into a format that SARGE can
% work with.

clear; close all
% Load .mat file into workspace
projRootPath = 'L:\My Drive\PROJECTS\STNstim_GPrecord\';
loadDirFullPath = 'D:\PROJECTS\STNstim_GPrecord\Data Acquisition\013020mat\';
matfn = 'jdsktestv2-200130-124237.mat'
streamStr = 'SUNx';

matFullPath = [loadDirFullPath matfn];
load(matFullPath);

% assume projRootPath was specified by "initProjResources.m"
% otherwise, specify it here:
% projRootPath = L:\My Drive\PROJECTS\STNstim_GPrecord\Data Processing\ArtifactSubtract;
saveDirFullPath = [projRootPath 'Data Processing\ArtifactSubtract\'];
% savePath = [saveDirFullPath matfn];



%% View file and get intervals to separate stim epocs in data (multiple DBS
% conditions)
streamData = data.streams.(streamStr).data(1,:);
fs = data.streams.(streamStr).fs;
t = (0:(size(streamData,2)-1)) * (1/fs);
figure; plot(t, streamData); grid minor
xlabel('time (seconds)')



%% Specify time points to break this data up if necessary
% default (i.e. no multiple blocks): [0, t(end)];

timePoints = [0, 70, 140, 200, 260, 320, t(end)]

% Now convert those times to samples
nBlocks = numel(timePoints) - 1;
blockSamps = zeros(nBlocks, 2);
for iBlk = 1:nBlocks
    blockSamps(iBlk,1) = round(timePoints(iBlk) * fs);
    blockSamps(iBlk,2) = round(timePoints(iBlk+1) * fs) - 1;
    
end
% Assure correct samples at beginning and end
blockSamps(1,1) = 1;
blockSamps(end,end) = length(streamData);



%% Parse matfile into seperate epoc blocks and channels, save them in
% current folder

saveParsedStream(matFullPath, streamStr, blockSamps, saveDirFullPath);
