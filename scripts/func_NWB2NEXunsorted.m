function func_NWB2NEXunsorted(nwb, spkChans, nexFullPath)
% script for generating a .nex file holding spike-filtered data ready for
% sorting in Plexon Offline sorter


% % clear; close all
% chanList = [64, 65, 66, 68, 70, 71, 74, 78];
% spkChanIdx = [4, 7]; % index for the above chans list

% % Read in nwb file
% nwbpn = 'C:\DATAtemp\STNstim_GPrecord\Data Acquisition\Nebula-210726\';
% nwbfn = 'AzNeuralverTwo-210726-105928';
% nwb = nwbRead([nwbpn nwbfn '.nwb']);

% iCh = 7;
[nexpn, nexfn, ~] = fileparts(nexFullPath);

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



%% Create NEX file with all data
nChans = length(spkChans);
for iCh = 1:nChans
    chanLabel = ['_ch' num2str(spkChans(iCh))];
    nexFile = nexCreateFileData(Fs);

    nexFile = nexAddContinuous(nexFile, 0, Fs, valuesMV(:,iCh), ...
        [nexfn chanLabel]);


    writeNexFile(nexFile, [nexpn '\' nexfn chanLabel '_unsorted' '.nex']);
    
    clear nexFile
    
end




