function func_combineNWBandSorted2NEXanalyze(nwb, spkChans, DBSelec, nwbfullpath)
%% Read in sorted NEX file, and add in DBS times as nex Events

% clear; close all
[fullpath, filename, ~] = fileparts(nwbfullpath);
nwbfullpath = [fullpath '\' filename];

% nwbpn = 'C:\DATAtemp\STNstim_GPrecord\Data Acquisition\Nebula-210726\';
% nwbfn = 'AzNeuralverTwo-210726-105928';
% nwb = nwbRead([nwbfullpath '.nwb']);

% nChans = length(spkChans);
for iCh = spkChans
    nexFile = readNexFile([nwbfullpath '_ch' num2str(iCh) '_SORTED.nex']);
    
    pm = nwb.processing.get('DBSartifact_detect-clean').nwbdatainterface;
    tDBS_virtPre = pm.get('tDBS_virtPre').timestamps.load;
            tDBS = pm.get('tDBS').timestamps.load;
    tDBS_virtPos = pm.get('tDBS_virtPos').timestamps.load;
    
    % for tDBS
    dbsType = 'tDBS';
    name = [DBSelec '_' dbsType];
    nexFile = nexAddEvent(nexFile, tDBS, name);

    % tDBS_virtPre
    dbsType = 'tVirtPre';
    name = [DBSelec '_' dbsType];
    nexFile = nexAddEvent(nexFile, tDBS_virtPre, name);

    % for tDBS
    dbsType = 'tVirtPos';
    name = [DBSelec '_' dbsType];  
    nexFile = nexAddEvent(nexFile, tDBS_virtPos, name);

    writeNexFile(nexFile, [nwbfullpath '_ch' num2str(iCh) '_' DBSelec '_ANALYZE.nex']);  
        
    clear nexFile

end



% % Read in nex file with sorted results
% nexpn = 'L:\My Drive\PROJECTS\STNstim_GPrecord\Data Processing\NeuroExplorer\Unsorted\';
% nexfn = [nwbfn '_ch74_SORTED'];
% nexFile = readNexFile([nwbfullpath '.nex']);

% DBSelec = 'E3';



% % for each dbs epoch, need to add in an event for DBS timings
% ep = nwb.intervals_epochs;
% start_time = ep.start_time.data.load;
% stop_time = ep.stop_time.data.load;
% nBlks = length(start_time);

% pm = nwb.processing.get('DBSartifact_detect-clean').nwbdatainterface;
% tDBS_virtPre = pm.get('tDBS_virtPre').timestamps.load;
%         tDBS = pm.get('tDBS').timestamps.load;
% tDBS_virtPos = pm.get('tDBS_virtPos').timestamps.load;
% 
% % for tDBS
% dbsType = 'tDBS';
% % for iBlk = 1:nBlks
%     name = [DBSelec '_' dbsType];
% %     timestamps = tDBS((tDBS >= start_time(iBlk)) & ...
% %                       (tDBS < stop_time(iBlk)));
%     
%     nexFile = nexAddEvent(nexFile, tDBS, name);
%     
% % end
% 
% % tDBS_virtPre
% dbsType = 'tVirtPre';
% % for iBlk = 1:nBlks
%     name = [DBSelec '_' dbsType];
% %     timestamps = tDBS_virtPre((tDBS_virtPre >= start_time(iBlk)) & ...
% %                               (tDBS_virtPre < stop_time(iBlk)));
%     
%     nexFile = nexAddEvent(nexFile, tDBS_virtPre, name);
%     
% % end
% 
% 
% % for tDBS
% dbsType = 'tVirtPos';
% % for iBlk = 1:nBlks
% name = [DBSelec '_' dbsType];  
% nexFile = nexAddEvent(nexFile, tDBS_virtPos, name);
%     
% % end
% 
% 
% 
% writeNexFile(nexFile, [nexpn nexfn(1:end-7) '_ANALYZE.nex']);


%% Create NEX file with all data
% nChans = length(spkChans);
% for iCh = 1:nChans
%     chanLabel = ['_ch' num2str(spkChans(iCh))];
%     nexFile = nexCreateFileData(Fs);
% 
%     nexFile = nexAddContinuous(nexFile, 0, Fs, valuesMV(:,iCh), ...
%         [nexfn chanLabel]);
% 
% 
%     writeNexFile(nexFile, [nexpn '\' nexfn chanLabel '_unsorted' '.nex']);
%     
%     clear nexFile
%     
% end


% INPUT:
%   nexFile - nex file data structure created in nexCreateFileData
%   timestamps - vector of event timestamps in seconds
%   name - event name



