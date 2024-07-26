% script to run data processing pipeline for STNstim_GPrecord project,
% specifically on Nebula's data. 

% This script is intended to run on one TDTblock file at a time, with the
% end result being one or more NEX files (for each spike-channel of iterenst
% specified. Each NEX file would have spike times, and DBS times (real and
% virtual), allowing for PSTH analysis in NeuroExplorer

% NOTE: based on prior script of same name (effectively v1), but with more
% modularity to accomodate processing of Azula's data as well as Nebula's





%% CODE
clear; close all

iTDTblock = 'jdsktestv2-200124-161037.mat';
projRootPn = 'C:\DATAtemp\STNstim_GPrecord\';
acqPn = 'Data Acquisition\';
prcPn = 'Data Processing\';

% Load data tables with metadata for each block (BlockData.xlsx) and each 
% session composed of multiple blocks (SessionData.xlsx)
blkTab = readtable([projRootPn prcPn 'BlockData_Azula.xlsx']);
blkTab2 = readtable([projRootPn prcPn 'BlockData_Azula_250uA.xlsx']);
sessTab = readtable([projRootPn prcPn 'SessionData_Azula.xlsx']);
dataTab = join(blkTab, blkTab2, 'Keys', 'TDTblock');
dataTab = join(dataTab, sessTab);

iBlkTab = dataTab(strcmp(dataTab.TDTblock, iTDTblock),:);

spkChans = str2num(iBlkTab.spkChans{:}); % Specifies which channels to process
DBSelec = iBlkTab.DBSelec{:};

% CONSTANTS
TDTblock = iBlkTab.subSelectName{:};
TDTtank = [iBlkTab.TDTtank{:} '\'];

cfg.artifact.threshold = -1e-4;
cfg.artifact.procSpkChans = spkChans;
cfg.artifact.artDetect_offest = -1; % samples
cfg.artifact.artDetect_blankWindow = [-4, 8]; % samples around DBS pulse
cfg.artifact.spkBandPass = [200, 8000]; % Hz
cfg.artifact.spkStore = 'SUNx';


%%  Convert TDTblock to NWB file, with all stream data stored under
% /acquisition section of NWB file
% function based on "script_TDTblock2NWB.m"

% streamParams(1).name = 'Wav2';
% streamParams(1).data_unit = 'Volts';
% streamParams(1).data_conversion = 1/(1e9);
% streamParams(2).name = 'pSUg';
% streamParams(2).data_unit = '<missing>';
% streamParams(2).data_conversion = 1;
% streamParams(3).name = 'Dbs1';
% streamParams(3).data_unit = '<missing>';
% streamParams(3).data_conversion = 1;
% streamParams(4).name = 'RSn1';
% streamParams(4).data_unit = '<missing>';
% streamParams(4).data_conversion = 1;

% nwb = util.TDTmat2NWB([projRootPn acqPn TDTtank TDTblock], streamParams);
nwb = util.TDTmat2NWB([projRootPn acqPn TDTtank TDTblock]);

% % In the case that Dbs1 was not used
% if ~any(strcmp(nwb.acquisition.keys, 'Dbs1'))
%     % First get "Dbs1" by stealing a channel from Wav2, and upsample it
%     sampFactor = 2;
%     wav2 = nwb.acquisition.get('Wav2');
%     chans = 68;
% 
%     % Get original timestamps
%     N = size(wav2.data, 1);
%     t0 = wav2.starting_time;
%     fs = wav2.starting_time_rate;
%     tst = t0 + ((0:(N-1))) * (1/fs);
% 
%     % Get Resampled timestamps
%     N_res = N * sampFactor;
%     t0_res = t0 / sampFactor;
%     fs_res = fs * sampFactor;
%     tst_res = ((1/fs_res) * (0:(N_res-1))) + t0_res;
% 
%     nChans = length(chans);
%     raw = zeros(N_res, nChans);
%     for iCh = 1:nChans
%         data = double(wav2.data(:,chans(iCh)));
%         data_res = interp1(tst, data, tst_res);
%         data_res(1) = data_res(2);
%         data_res(end) = data_res(end-1); % correct for any NaNs output by interp1...
%         raw(:,iCh) = data_res;
% 
%     end




    

%     data = TDTbin2mat(blockpath, 'STORE', tdtstoreString);

    % Then save it in acquisition segment of nwb file
    
   
% Choose a particular timeseries in acquisition for DBS detection, saved as
% a n x 1 vector voltage trace



        % Store data as an ElectricalSeries object
            tdtstoreString = 'Dbs1';
            
        tSeries = nwb.acquisition.get(iBlkTab.stimDetectStore);    
        dataDbs = tSeries.data(:,iBlkTab.stimDetectChan);

        SeriesName = tdtstoreString;
%         time_series = types.core.TimeSeries( ...
%             'data', dataDbs, ...
%             'data_resolution', -1.0, ... % default for unknown
%             'data_continuity', 'continuous', ...
%             'data_unit', '<missing>', ...
%             'data_conversion', 1, ...
%             'starting_time', 0, ... % seconds
%             'starting_time_rate', 24414.0625);
        time_series = types.core.TimeSeries( ...
            'data', dataDbs, ...
            'data_resolution', tSeries.data_resolution, ... % default for unknown
            'data_continuity', tSeries.data_continuity, ...
            'data_unit', tSeries.data_unit, ...
            'data_conversion', tSeries.data_conversion, ...
            'starting_time', tSeries.starting_time, ... % seconds
            'starting_time_rate', tSeries.starting_time_rate);

        nwb.acquisition.set(SeriesName, time_series);

%     end
    
% end

% write and read back an NWB file
% write file
nwbfn = TDTblock; % gets rid of the .mat part of the string
nwbfn = util.remextension(nwbfn); % custom code to remove any extension chars
nwbfullpath = [projRootPn prcPn TDTtank nwbfn '.nwb'];
disp(['Exporting ' nwbfn '.nwb...']);
tic
nwbExport(nwb, nwbfullpath)
toc
disp('DONE!')

% Save RAM memory by reloading the nwb file with lazy loading
clear nwb
nwb = nwbRead(nwbfullpath);

% test read
% nwbTest = nwbRead(nwbfullpath, 'ignoreCache');



%% Perform DBS-artifact detection and cleaning, with the following
% stored under /processing of NWB file: DBS times(real and virtual), and
% spkClean continuous data
% fucntion based on "script_removeArtifacts_NWB.m"

nwb = func_removeArtifacts_NWB(nwb, cfg.artifact);

% Save updates to NWB file
nwb = util.nwbOverwrite(nwb, nwbfullpath);
disp([nwbfn '.nwb updated with processing for DBS artifact detecting and cleaning.']);



%% Generate NEX file with unsorted continuous data for each specified
% channel of interest

func_NWB2NEXunsorted(nwb, spkChans, nwbfullpath);
% After sorting is done on each NEX-unsorted file in Plexon Offline sorter,
% Add spike data to NWB file

% Based on spike data and DBS time data, generate NEX file for analysis on
% NeuroExplorer for each channel of interest



%% Generate finished NEX file with sorted spike times and DBS event times

% NOTE!!! Put a break at this line so that you have a chance to sort any
% files in PLEXON first!!

disp('Perform all Spike Sorting in PLEXON now!!')
pause
func_combineNWBandSorted2NEXanalyze(nwb, spkChans, DBSelec, nwbfullpath)
disp('NEX files ready for analysis!')



