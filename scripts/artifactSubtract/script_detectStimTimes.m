% script to detect DBS stim times based on the "Dbs1" store in the TDT
% block recorded, now in nwb file format

clear; clc;
% Load data
nwbpn = 'C:\DATAtemp\STNstim_GPrecord\Data Acquisition\Nebula-210726\';
nwbfn = 'AzNeuralverTwo-210726-105657';
nwb = nwbRead([nwbpn nwbfn '.nwb']);

Dbs1 = nwb.acquisition.get('Dbs1');
dbsTrace = double(Dbs1.data.load);
fs = Dbs1.starting_time_rate;
tStamp = (1/fs) * (0:(length(dbsTrace) - 1));
figure; plot(tStamp, dbsTrace)


% Detect stim times

% 
% clear
% % Load TDT file into workspace
% % loadDirFullPath = 'L:\My Drive\PROJECTS\STNstim_GPrecord\Data Acquisition\Sokka\';
% loadDirFullPath = 'C:\DATAtemp\STNstim_GPrecord\Data Processing\Artifact Subtraction\';
% 
% matfn = 'AzNeuralverTwo-210726-105657_Dbs1.mat';
% [~, matfn, ~] = fileparts(matfn); % remove any extensions 
% 
% % streamStr = 'SUNx';
% 
% % matFullPath = [loadDirFullPath TDTblock];
% % load(matFullPath);8
% 
% % detectionChan = 'SUNx_ch3';
% disp('Loading mat file...')
% load([loadDirFullPath matfn]);
% disp('Done loading mat file')
% 
% 
% % assume projRootPath was specified by "initProjResources.m"
% % otherwise, specify it here:
% % projRootPath = 'L:\My Drive\PROJECTS\STNstim_GPrecord\';
% saveDirFullPath = loadDirFullPath;
% % savePath = [saveDirFullPath matfn];
% 
% 
% % CONSTANTS
% 
% FS = 24414.0625;
% % BLK = 1;


% 
% %%
% % streamStr = 'SUNx';
% close all
% streamStr = 'SUNx';
% 
% % recElec = 1; % manually change this based on observation above code
% 
% eval(['dbsTrace = ' 'Dbs1' ';'])
% 
% fs = FS;
% tStamp = (1/fs) * (0:(length(dbsTrace) - 1));
% 
% figure; plot(tStamp, dbsTrace)
% 


%%

THRESH = -0.01; % manually change this basded on observation
x = dbsTrace; 

if THRESH < 0 % if thresh is a negative crossing, flip the data over 
    x = -x; 
    THRESH = -THRESH;
    
end

% Essentially need to pick out values which exceed threshold with the condition that the previous value  
% needs to be below the threshold
idxl = x >= THRESH;
idxl(1) = 0;
idx = find(idxl);
yest = x(idx-1) < THRESH; 
idxStim = idx(yest); % Final output
idxStim = idxStim - 1; % correct by shifting stim onset one sample back
stimTime = tStamp(idxStim);

disp(['detected stims: ' num2str(length(stimTime))]);
% should be just one bin; more bins means that multiple parts of the pulse
% are getting detected 
figure; histogram(diff(stimTime));



%% Update the current nwb file to include DBS stim times

stim_series_DBS = types.core.TimeSeries(...
    'data', ones(length(stimTime), 1), ...
    'data_unit', 'N/A', ...
    'description', 'DBS event timings found in timestamps field, the data field is just placeholder values', ...
    'timestamps', stimTime);

ecephys_module = types.core.ProcessingModule(...
    'description', 'holds extracellular electrophysiology data');
ecephys_module.nwbdatainterface.set('LFP', ...
    types.core.LFP('lfp', electrical_series));
nwb.processing.set('ecephys', ecephys_module);




nwb.processing.set('DBStimes', stim_series_DBS);

% % This fig takes a while to plot... but useful: plots dotted vertical lines
% % wherever pulses detected
% figure; ax = axes; plot(tStamp, dbsTrace); hold on
% for i = 1:numel(idxStim)
% %     plot([idxStim(i), idxStim(i)], ax.YLim, '--r')
%     plot([stimTime(i), stimTime(i)], ax.YLim, '--r')
% 
%     
% end



% Parse matfile into seperate epoc blocks and channels, save them in
% current folder

% Also save the detected stimTimes
% Note: This replaces the detection step that I had originally done in
% SARGE for Azula's data...
disp('Saving stimtimes...')
saveFullPathStimtimes = [saveDirFullPath matfn '_stimtimes.mat'];
save(saveFullPathStimtimes, 'stimTime');   
disp('Done saving stimtimes!')


