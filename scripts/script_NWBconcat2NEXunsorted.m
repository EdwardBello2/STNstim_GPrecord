% script for generating a .nex file holding spike-filtered data ready for
% sorting in Plexon Offline sorter


clear; close all
chanList = [64, 65, 66, 68, 70, 71, 74, 78];
spkChanIdx = [4, 7]; % index for the above chans list

% Read in nwb file
nwbpn = 'C:\DATAtemp\STNstim_GPrecord\Data Acquisition\Nebula-210726\';
nwbfn = 'StnDbsNeb_210726-105657';
nwb = nwbRead([nwbpn nwbfn '.nwb'])

% iCh = 7;

% Extract necessary data
spkClean = nwb.processing.get('DBSartifact_detect-clean').nwbdatainterface.get('spkClean');
Fs = spkClean.starting_time_rate; % may need to derive Fs differently later..
t0 = spkClean.starting_time;
values = spkClean.data.load;
values = values * spkClean.data_conversion; % Get data in volts
valuesMV = values * 1000; % Convert from volts to millivolts as NEX requires

% A little trick to curtail writeNexFile's habit of truncating upper and
% lower voltages on spike waveforms (perhaps a data compression measure?)
valuesMV(1,:) = 2 * max(max(valuesMV));

%% Get sample info of each data fragment

tst = spkClean.timestamps.load;
epochIntervals = nwb.intervals_epochs;
idxStart = epochIntervals.vectordata.get('spkClean_idx_start').data.load;
 idxStop = epochIntervals.vectordata.get('spkClean_idx_stop').data.load;
tstStart = tst(idxStart);



% %% Create stitched-together version of session data with zeros filling in intervening time
% % Iteratively append zeros followed by real data
% 
% Ch = 7;
% dataMV = valuesMV(:,Ch);
% 
% nBlocks = length(idxStart);
% for iBlk = 1:nBlocks
%     itst = tst(idxStart(iBlk):idxStop(iBlk));
%     idataMV = dataMV(idxStart(iBlk):idxStop(iBlk));
%     % fill timestamps
%     
%     
% end
% 
% 

%% Create NEX file with all data

for iCh = spkChanIdx
    chanLabel = num2str(chanList(iCh));
    nexFile = nexCreateFileData(Fs);

    % specify start time (t(1)), digitizing frequency (Fs), data (x2) and name
    % for iCh = 1:nChans
    %     nexFile = nexAddContinuous(nexFile, t0, Fs, valuesMV(:,iCh), ['ch' num2str(iCh)]);

    % end
    nexFile.contvars{1}.name = 'test';
    nexFile.contvars{1}.ADFrequency = Fs;
    nexFile.contvars{1}.timestamps = tstStart;
    nexFile.contvars{1}.fragmentStarts = idxStart;
    nexFile.contvars{1}.data = valuesMV(:,iCh);
    nexFile.contvars{1}.varVersion = 100;
    nexFile.tbeg = tst(1);
    nexFile.tend = tst(end);




    % save nex file (assuming nex directory in Windows 7)
    nexpn = 'L:\My Drive\PROJECTS\STNstim_GPrecord\Data Processing\NeuroExplorer\Unsorted\';
    nexfn = nwbfn;

    writeNexFile(nexFile, [nexpn nexfn '_ch' chanLabel '_unsorted' '.nex']);
    
    clear nexFile
    
end

% 
% 
% %% Read in sorted NEX file, and add in DBS times as nex Events
% 
% % Read in nex file with sorted results
% nexpn = 'L:\My Drive\PROJECTS\STNstim_GPrecord\Data Processing\NeuroExplorer\Unsorted\';
% nexfn = 'testSess7_sorted';
% nexFile = readNexFile([nexpn nexfn '.nex']);
% 
% % Add as event in nex file
% TS = nwb.processing.get('DBSartifact_detect-clean').nwbdatainterface.get('tDBS');
% timestamps = TS.timestamps.load;
% nexFile = nexAddEvent(nexFile, timestamps, 'DBS');
% 
% 
% writeNexFile(nexFile, [nexpn nexfn num2str(iCh) '_analyze.nex']);
% 
% 
% 
% 
% 
% % INPUT:
% %   nexFile - nex file data structure created in nexCreateFileData
% %   timestamps - vector of event timestamps in seconds
% %   name - event name
% 


