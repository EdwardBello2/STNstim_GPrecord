%% Read in sorted NEX file, and add in DBS times as nex Events

clear; close all

nwbpn = 'C:\DATAtemp\STNstim_GPrecord\Data Acquisition\Nebula-210726\';
nwbfn = 'AzNeuralverTwo-210726-105928';
nwb = nwbRead([nwbpn nwbfn '.nwb']);


% Read in nex file with sorted results
nexpn = 'L:\My Drive\PROJECTS\STNstim_GPrecord\Data Processing\NeuroExplorer\Unsorted\';
nexfn = [nwbfn '_ch74_SORTED'];
nexFile = readNexFile([nexpn nexfn '.nex']);

dbsEpochLabel = 'E3';



% % for each dbs epoch, need to add in an event for DBS timings
% ep = nwb.intervals_epochs;
% start_time = ep.start_time.data.load;
% stop_time = ep.stop_time.data.load;
% nBlks = length(start_time);

pm = nwb.processing.get('DBSartifact_detect-clean').nwbdatainterface;
tDBS_virtPre = pm.get('tDBS_virtPre').timestamps.load;
        tDBS = pm.get('tDBS').timestamps.load;
tDBS_virtPos = pm.get('tDBS_virtPos').timestamps.load;

% for tDBS
dbsType = 'tDBS';
% for iBlk = 1:nBlks
    name = [dbsEpochLabel '_' dbsType];
%     timestamps = tDBS((tDBS >= start_time(iBlk)) & ...
%                       (tDBS < stop_time(iBlk)));
    
    nexFile = nexAddEvent(nexFile, tDBS, name);
    
% end

% tDBS_virtPre
dbsType = 'tVirtPre';
% for iBlk = 1:nBlks
    name = [dbsEpochLabel '_' dbsType];
%     timestamps = tDBS_virtPre((tDBS_virtPre >= start_time(iBlk)) & ...
%                               (tDBS_virtPre < stop_time(iBlk)));
    
    nexFile = nexAddEvent(nexFile, tDBS_virtPre, name);
    
% end


% for tDBS
dbsType = 'tVirtPos';
% for iBlk = 1:nBlks
    name = [dbsEpochLabel '_' dbsType];
%     timestamps = tDBS_virtPos((tDBS_virtPos >= start_time(iBlk)) & ...
%                               (tDBS_virtPos < stop_time(iBlk)));
    
    nexFile = nexAddEvent(nexFile, tDBS_virtPos, name);
    
% end



writeNexFile(nexFile, [nexpn nexfn(1:end-7) '_ANALYZE.nex']);





% INPUT:
%   nexFile - nex file data structure created in nexCreateFileData
%   timestamps - vector of event timestamps in seconds
%   name - event name



