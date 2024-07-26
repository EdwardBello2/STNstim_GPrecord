% script for generating a .nex file holding spike-filtered data ready for
% sorting in Plexon Offline sorter


clear; close all

% Read in nwb file
nwbpn = 'C:\DATAtemp\STNstim_GPrecord\Data Acquisition\Nebula-210726\';
nwbfn = 'AzNeuralverTwo-210726-110153';
nwb = nwbRead([nwbpn nwbfn '.nwb'])

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

% Create NEX file with all data
nexFile = nexCreateFileData(Fs);

% specify start time (t(1)), digitizing frequency (Fs), data (x2) and name
nChans = size(values, 2);
for iCh = 1:nChans
    nexFile = nexAddContinuous(nexFile, t0, Fs, valuesMV(:,iCh), ['ch' num2str(iCh)]);

end


% save nex file (assuming nex directory in Windows 7)
nexpn = 'L:\My Drive\PROJECTS\STNstim_GPrecord\Data Processing\NeuroExplorer\Unsorted\';
nexfn = nwbfn;

writeNexFile(nexFile, [nexpn nexfn '.nex']);

