% script for saving SARGE-compatible data from NWBsession data, one block
% at a time...



%% Read in NWB file for whole recording session
pn = 'C:\DATAtemp\STNstim_GPrecord\Data Acquisition\Nebula-210726\';
nwbfn = 'AzNeuralverTwo-210726-105657';
nwb = nwbRead([pn nwbfn '.nwb'],'ignoreCache');


%% Read in one block of data from the session (pre DBS, DBS on, post DBS)

% Get info from epochs into a table format
streamStr = 'Dbs1';
tseries = nwb.acquisition.get(streamStr)

Dbs1 = tseries.data.load;


% Save as .mat file in proper format
savepn = 'C:\DATAtemp\STNstim_GPrecord\Data Processing\Artifact Subtraction\';
save([savepn nwbfn '_' streamStr], 'Dbs1');