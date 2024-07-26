% test script to derive lfp data from an epoch around voltage trace

clear
%% CONSTANTS

SESSION_ID = 'TremoLfpDBS-190927-100155';

% bipolar electrode pair:
ePair = [1, 3];
elecLabels = {'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'};

nwbProcSuffix = ['LFP_' elecLabels{ePair(1)} '-' elecLabels{ePair(2)} '_washin_exp'];
nwbReadAcquPn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';
nwbWritePn = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Processing\NWBdata\';

tic
%%  Get bipolar lfp channel C0-C2 data for washin
clear raw epochIntervals

% read MASTER nwb file (with raw-acquisition data) 
nwbAcqu = nwbRead([nwbReadAcquPn SESSION_ID '_MASTER.nwb']);
raw = nwbAcqu.acquisition.get('rawVoltageTrace').deref();

% Get relevant interval times for subselection of data
epochIntervals = nwbAcqu.intervals_epochs.deref();
eT = nwbnrtl.util.readEpochsTable(epochIntervals);
epoch_tag = 'harWashinExp';
epochRow = strcmp(eT.tag, epoch_tag);
interval = [eT{epochRow, 'start_time'}, eT{epochRow, 'stop_time'}];

% Get subselected lfp data and timestamps
data = nwbnrtl.util.loadTimeSeriesData(raw, interval, [], [1, 3]);    
tstamps = nwbnrtl.util.loadTimeSeriesTimestamps(raw, interval);

dataDiff = diff(data, 1, 2);
fs = raw.starting_time_rate;



%% LFP-filter and downsample to ~ 2000 samples/sec

downsampFactor = 24; % get it from 48kHz to ~ 2kHz

fc = [0.1, 500];

% Remove 60 Hz line noise
% dataDiff = filtLineNoise(dataDiff, fs, 5);

[b, a] = butter(2, (fc / (fs/2)), 'bandpass');
dataFilt = filtfilt(b, a, dataDiff);
dataLfp = downsample(dataFilt, downsampFactor);
fsLfp = fs/downsampFactor;
tstampsLfp = downsample(tstamps, downsampFactor);

% 
% 
% % run spectrogram on it
% window = 2^10;
% noverlap = floor(0.25*window);
% f1 = figure; 
% spectrogram(dataLfp, window, noverlap, [], fsLfp, 'yaxis');
% ax1 = gca;
% ax1.CLimMode = 'manual';
% ax1.CLim = [-130 -110]; % [-156.5356 -77.2316]
% ax1.YLim = [-0.4967 100]; % [-0.4967 509.1230]
% f1.Position = [2599 475 1379 420];
% 
% % Look at PSD's over time
% t = (1/fsLfp) * (0:(length(dataLfp)-1));
% secWindow = floor(300 * fsLfp)
% % plot psd from 1 to 300 seconds
% % [pxx, f] = pwelch(dataLfp(1:secWindow), fsLfp);
% f2 = figure;
% ax2 = axes;
% % pwelch(dataLfp(1:secWindow), fsLfp);
% pwin = 2^11;
% poverlap = 2^5;
% pnfft = pwin;
% [pxx, f] = pwelch(dataLfp(1:secWindow), pwin, poverlap, pnfft, fsLfp);
% plot(f, 10*log10(pxx))
% ylabel('PSD (dB/Hz)')
% xlabel('Frequency (Hz)')
% ax2.XLim = [0, 55];



%% generate lfp nwb file

id = genUUID;

% Init NWB
nwbProc = NwbFile(...
    'session_start_time', nwbAcqu.session_start_time, ...
    'identifier', id, ...
    'session_description', 'missing', ...
    'timestamps_reference_time', nwbAcqu.session_start_time);



% Fill "general" section of nwb file, constants set inside sub-function 
nwbProc = nwbfill_general(nwbProc);
nwbProc.general_session_id = SESSION_ID;
[~,fname,~] = fileparts(mfilename('fullpath'));
nwbProc.general_source_script_file_name = fname;

% Fill "electrode" information of nwb file, constants set inside sub-function
nwbProc = nwbfill_electrode(nwbProc);


% specify the electrode(s) that this signal came from according to nwb
electrodes_object_view = types.untyped.ObjectView( ...
    '/general/extracellular_ephys/electrodes');

electrode_table_region = types.hdmf_common.DynamicTableRegion( ...
    'table', electrodes_object_view, ...
    'description', 'reference to the two channels from which bipolar LFP was derived.', ...
    'data', [0, 2]); % nwb uses zero as first index...


% store lfps as electrical series object in nwb under Processing
electrical_series = types.core.ElectricalSeries( ...
    'comments', 'Bipolar LFP derived from "rawVoltageTrace" in acquisision, using "script_genNWB_proc_LFP.m"', ...
    'data', dataLfp, ...
    'data_resolution', -1.0, ... % default for unknown
    'data_unit', 'Volts', ...
    'description', 'data stored as single-precision matlab array', ...
    'electrodes', electrode_table_region, ... % points to value in electrode table
    'starting_time', tstampsLfp(1), ... % seconds
    'starting_time_rate', fsLfp, ...
    'timestamps', tstampsLfp);
    
    

ecephys_module = types.core.ProcessingModule(...
    'description', 'holds extracellular electrophysiology data');
ecephys_module.nwbdatainterface.set('LFP', ...
    types.core.LFP('lfp', electrical_series));
nwbProc.processing.set('ecephys', ecephys_module);



%% write nwb file for lfp data

nwbExport(nwbProc, [nwbWritePn SESSION_ID '_' nwbProcSuffix '.nwb'])
% nwbExport(nwb, [NWBPATH 'TEST' '.nwb'])

disp('DONE EXPORTING NWB!');
toc



%% SUB-FUNCTIONS

function nwb = nwbfill_general(nwb)

nwb.general_data_collection = {'test'};
nwb.general_experiment_description = {'test'};
nwb.general_experimenter = {'Bello', 'Blumenfeld'};
% nwb.general_extracellular_ephys = 'test' % need to fix
% nwb.general_extracellular_ephys_electrodes = 'test' % need to fix 
nwb.general_institution = {'University of Minnesota'};
nwb.general_keywords = {'harmaline', 'touchscreen', 'lfp'};
nwb.general_lab = {'NRTL'};
nwb.general_notes = {'missing'};
nwb.general_pharmacology = {'missing'};
nwb.general_protocol = {'missing'};
nwb.general_subject = types.core.Subject(...
    'date_of_birth', '2001-01-01 00:00:00', 'description', 'nickname: Uva, source: breeder', ...
    'sex', 'F', 'species', 'macacca mulatta', 'subject_id', '13LP1');

nwb.general_devices.set('PZ5', types.core.Device());

end

function nwb = nwbfill_electrode(nwb)

deviceName = '2017_Numed-002';
nwb.general_devices.set(deviceName, types.core.Device());

% specify subgrouping of electrodes in recording device (here there's just
% one group, that is, all the DBS electrodes)
elecGroupName = 'DBSelectrodes';
nwb.general_extracellular_ephys.set(elecGroupName, ...
    types.core.ElectrodeGroup(...
    'description', 'linear array of annular electrodes', ...
    'location', 'lateral thalamus', ...
    'device', types.untyped.SoftLink(['/general/devices/', deviceName])));


% % set up electrode information table
% group_object_view = types.untyped.ObjectView( ...
%     ['/general/extracellular_ephys/' elecGroupName]);

% tableKeys = {'label', 'impedance', 'electrodeGroup'};
tableKeys = {'label', 'impedance'};
labels = {'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'}'; 
% nElecs = numel(labels);
% id = [1:nElecs]';
imp = NaN(8,1);
tab = table(labels, imp);
% tab = table(labels(1), {NaN}, group_object_view);
% tab = table(1, labels(1), {NaN});
% 
% for iE = 2:nElecs
% %     tab = [tab; table(labels(iE), {NaN}, group_object_view)];
%         tab = [tab; table(iE, labels(iE), {NaN})];
% 
% end
% tab.Properties.DimensionNames = {'electrode', 'property'};
tab.Properties.VariableNames = tableKeys;
% tab.Properties.VariableDescriptions = {'naming convenction', 'abs imp @ 1kHz'};
tab.Properties.Description = 'information for use within nwb file';

elecDynTable = util.table2nwb(tab);
nwb.general_extracellular_ephys_electrodes = elecDynTable;

end

function id = genUUID
temp = java.util.UUID.randomUUID;
myuuid = temp.toString;
id = char(myuuid);
% fprintf('%s\n', char(myuuid))

end

