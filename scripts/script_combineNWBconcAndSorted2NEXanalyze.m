% % script for generating a .nex file holding spike-filtered data ready for
% % sorting in Plexon Offline sorter
% 
% 
% clear; close all
% chansList = [64, 65, 66, 68, 70, 71, 74, 78];
% spkChanIdx = [4, 7]; % index for the above chans list
% 
% % Read in nwb file
% nwbpn = 'C:\DATAtemp\STNstim_GPrecord\Data Acquisition\Nebula-210726\';
% nwbfn = 'StnDbsNeb_210726-105657';
% nwb = nwbRead([nwbpn nwbfn '.nwb'])
% 
% % iCh = 7;
% 
% % Extract necessary data
% spkClean = nwb.processing.get('DBSartifact_detect-clean').nwbdatainterface.get('spkClean');
% Fs = spkClean.starting_time_rate; % may need to derive Fs differently later..
% t0 = spkClean.starting_time;
% values = spkClean.data.load;
% values = values * spkClean.data_conversion; % Get data in volts
% valuesMV = values * 1000; % Convert from volts to millivolts as NEX requires
% 
% % A little trick to curtail writeNexFile's habit of truncating upper and
% % lower voltages on spike waveforms (perhaps a data compression measure?)
% valuesMV(1,:) = 2 * max(max(valuesMV));
% 
% %% Get sample info of each data fragment
% 
% tst = spkClean.timestamps.load;
% epochIntervals = nwb.intervals_epochs;
% idxStart = epochIntervals.vectordata.get('spkClean_idx_start').data.load;
%  idxStop = epochIntervals.vectordata.get('spkClean_idx_stop').data.load;
% tstStart = tst(idxStart);
% 
% 
% 
% % %% Create stitched-together version of session data with zeros filling in intervening time
% % % Iteratively append zeros followed by real data
% % 
% % Ch = 7;
% % dataMV = valuesMV(:,Ch);
% % 
% % nBlocks = length(idxStart);
% % for iBlk = 1:nBlocks
% %     itst = tst(idxStart(iBlk):idxStop(iBlk));
% %     idataMV = dataMV(idxStart(iBlk):idxStop(iBlk));
% %     % fill timestamps
% %     
% %     
% % end
% % 
% % 
% 
% %% Create NEX file with all data
% 
% for iCh = spkChanIdx
%     chanLabel = num2str(chanList(iCh));
%     nexFile = nexCreateFileData(Fs);
% 
%     % specify start time (t(1)), digitizing frequency (Fs), data (x2) and name
%     % for iCh = 1:nChans
%     %     nexFile = nexAddContinuous(nexFile, t0, Fs, valuesMV(:,iCh), ['ch' num2str(iCh)]);
% 
%     % end
%     nexFile.contvars{1}.name = 'test';
%     nexFile.contvars{1}.ADFrequency = Fs;
%     nexFile.contvars{1}.timestamps = tstStart;
%     nexFile.contvars{1}.fragmentStarts = idxStart;
%     nexFile.contvars{1}.data = valuesMV(:,iCh);
%     nexFile.contvars{1}.varVersion = 100;
%     nexFile.tbeg = tst(1);
%     nexFile.tend = tst(end);
% 
% 
% 
% 
%     % save nex file (assuming nex directory in Windows 7)
%     nexpn = 'L:\My Drive\PROJECTS\STNstim_GPrecord\Data Processing\NeuroExplorer\Unsorted\';
%     nexfn = nwbfn;
% 
%     writeNexFile(nexFile, [nexpn nexfn '_ch' chanLabel '_unsorted' '.nex']);
%     
%     clean nexFile
%     
% end



%% Read in sorted NEX file, and add in DBS times as nex Events

clear; close all
% Read in nex file with sorted results
nexpn = 'L:\My Drive\PROJECTS\STNstim_GPrecord\Data Processing\NeuroExplorer\Unsorted\';
nexfn = 'StnDbsNeb_210726-105657_ch74_SORTED';
nexFile = readNexFile([nexpn nexfn '.nex']);

dbsEpochLabel = {'E5', 'E3', 'E1', 'E6', 'E4', 'E2', ...
    'E3-E1', 'E6-E4', 'E2-E3', 'E1-E2', 'E4-E5', 'E5-E6'};

nwbpn = 'C:\DATAtemp\STNstim_GPrecord\Data Acquisition\Nebula-210726\';
nwbfn = 'StnDbsNeb_210726-105657';
nwb = nwbRead([nwbpn nwbfn '.nwb']);

% for each dbs epoch, need to add in an event for DBS timings
ep = nwb.intervals_epochs;
start_time = ep.start_time.data.load;
stop_time = ep.stop_time.data.load;
nBlks = length(start_time);

pm = nwb.processing.get('DBSartifact_detect-clean').nwbdatainterface;
tDBS_virtPre = pm.get('tDBS_virtPre').timestamps.load;
        tDBS = pm.get('tDBS').timestamps.load;
tDBS_virtPos = pm.get('tDBS_virtPos').timestamps.load;

% for tDBS
dbsType = 'tDBS';
for iBlk = 1:nBlks
    name = [dbsEpochLabel{iBlk} '_' dbsType];
    timestamps = tDBS((tDBS >= start_time(iBlk)) & ...
                      (tDBS < stop_time(iBlk)));
    
    nexFile = nexAddEvent(nexFile, timestamps, name);
    
end

% tDBS_virtPre
dbsType = 'tVirtPre';
for iBlk = 1:nBlks
    name = [dbsEpochLabel{iBlk} '_' dbsType];
    timestamps = tDBS_virtPre((tDBS_virtPre >= start_time(iBlk)) & ...
                              (tDBS_virtPre < stop_time(iBlk)));
    
    nexFile = nexAddEvent(nexFile, timestamps, name);
    
end


% for tDBS
dbsType = 'tVirtPos';
for iBlk = 1:nBlks
    name = [dbsEpochLabel{iBlk} '_' dbsType];
    timestamps = tDBS_virtPos((tDBS_virtPos >= start_time(iBlk)) & ...
                              (tDBS_virtPos < stop_time(iBlk)));
    
    nexFile = nexAddEvent(nexFile, timestamps, name);
    
end



writeNexFile(nexFile, [nexpn nexfn(1:end-7) '_ANALYZE.nex']);





% INPUT:
%   nexFile - nex file data structure created in nexCreateFileData
%   timestamps - vector of event timestamps in seconds
%   name - event name



