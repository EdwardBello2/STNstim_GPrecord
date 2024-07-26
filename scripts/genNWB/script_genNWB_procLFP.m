% script for generating whole session long LFP, saving as nwb processed
% file (to be referenced in the MASTER nwb file)


% Load raw ephys NWB
%% CONSTANTS

% SESSION_ID = 'TremoLfpDBS-190927-100155';

% % bipolar electrode pair:
% ePair = [5, 7];
% elecLabels = {'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'};
SESSION_ID = 'TremoLfpDBS-191115-100127';
nwbProcSuffix = ['_proc_lfp'];
nwbReadAcquPn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';
nwbWritePn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';

tic
%%  Get bipolar lfp channel C0-C2 data for washin
clear raw epochIntervals

% read MASTER nwb file (with raw-acquisition data) 
nwbAcqu = nwbRead([nwbReadAcquPn SESSION_ID '_MASTER.nwb']);
raw = nwbAcqu.acquisition.get('rawVoltageTrace').deref();
fsRaw = raw.starting_time_rate;

decFactor = 24;
fsLfp = fsRaw / decFactor;


% Collect sample lims for the multiple tdt "files" represented in NWB
% all values to be found in the string within raw.description
sampLims = [1, 418291712; ...
            418291713, 506912768; ...
            506912769, 592953344; ...
            592953345, 720818176; ...
            720818177, 798003200; ...
            798003201, 883650560];
        
        
        
% Loop thru decimating raw ephys and filtering to LFP, storing in LFP
% matrix
lfpMat = [];
tstMat = [];
nFiles = length(sampLims);
nChans = 8;
tic
for iFile = 1:nFiles
    
    % Decimate and concatenate lfp data
    lfpCell = cell(1, nChans);
    for iCh = 1:nChans
        % load one channel of one file at a time
        intSamp = sampLims(iFile,:);
        data = nwbnrtl.util.loadTimeSeriesData(raw, ...
            [], [], iCh, ...
            'intervalSamples', sampLims(iFile,:));
        
        
        % decimate to new LFP sampling rate
        dataDec = decimate(double(data), decFactor, 'fir');
        lfpCell{iCh} = dataDec;  
        
        disp(['finished adding chan ' num2str(iCh) ' of file ' num2str(iFile)])
        
    end
    temp = [lfpCell{1}, lfpCell{2}, lfpCell{3}, lfpCell{4}, ...
    lfpCell{5}, lfpCell{6}, lfpCell{7}, lfpCell{8}];
    lfpMat = [lfpMat; temp];
    
    
    % Downsample and concatenate timestamps to match the above lfp
    tstRaw = nwbnrtl.util.loadTimeSeriesTimestamps(raw, ...
            'intervalSamples', sampLims(iFile,:));
    tstDec = downsample(tstRaw, decFactor);
    tstMat = [tstMat; tstDec];
    
    
end
clear data
toc


% % Get relevant interval times for subselection of data
% epochIntervals = nwbAcqu.intervals_epochs.deref();
% eT = nwbnrtl.util.readEpochsTable(epochIntervals);
% epoch_tag = 'harWashinExp';
% epochRow = strcmp(eT.tag, epoch_tag);
% interval = [eT{epochRow, 'start_time'}, eT{epochRow, 'stop_time'}];
% 
% % Get subselected lfp data and timestamps
% data = nwbnrtl.util.loadTimeSeriesData(raw, interval, [], [ePair(1), ePair(2)]);    
% tstamps = nwbnrtl.util.loadTimeSeriesTimestamps(raw, interval);
% 
% dataDiff = diff(data, 1, 2);
% fs = raw.starting_time_rate;
% 


% %% LFP-filter and downsample to ~ 2000 samples/sec
% 
% downsampFactor = 24; % get it from 48kHz to ~ 2kHz
% 
% fc = [0.1, 500];
% 
% % Remove 60 Hz line noise
% % dataDiff = filtLineNoise(dataDiff, fs, 5);
% 
% [b, a] = butter(2, (fc / (fs/2)), 'bandpass');
% dataFilt = filtfilt(b, a, dataDiff);
% dataLfp = downsample(dataFilt, downsampFactor);
% fsLfp = fs/downsampFactor;
% tstampsLfp = downsample(tstamps, downsampFactor);
% 
% 

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
    'data', [0:7]); % nwb uses zero as first index...


% store lfps as electrical series object in nwb under Processing
electrical_series = types.core.ElectricalSeries( ...
    'comments', 'monopolar LFP derived from "rawVoltageTrace" in acquisision, using "script_genNWB_procLFP.m"', ...
    'data', lfpMat, ...
    'data_resolution', -1.0, ... % default for unknown
    'data_unit', 'Volts', ...
    'description', 'data stored as single-precision matlab array', ...
    'electrodes', electrode_table_region, ... % points to value in electrode table
    'starting_time', tstMat(1), ... % seconds
    'starting_time_rate', fsLfp, ...
    'timestamps', tstMat);
    
    

ecephys_module = types.core.ProcessingModule(...
    'description', 'holds extracellular electrophysiology data');
ecephys_module.nwbdatainterface.set('LFP', ...
    types.core.LFP('lfp', electrical_series));
nwbProc.processing.set('ecephys', ecephys_module);



%% write nwb file for lfp data

nwbExport(nwbProc, [nwbWritePn SESSION_ID nwbProcSuffix '.nwb'])
% nwbExport(nwb, [NWBPATH 'TEST' '.nwb'])

disp('DONE EXPORTING NWB!');
toc

% end

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