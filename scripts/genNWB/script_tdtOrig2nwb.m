% test script for converting raw Acquisition data from tdt (and IMU &
% potentially camera data) into one consolidated NWB file
clear

%% Get the number of samples from each of the TDT blocks to be combined

SESSION_ID = 'testSession5';

tdtDataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';
TDTtank = 'Uva-191115';
IMUdata = 'TremoLfpDBS-191115-100127_IMU.csv';


% Specify the TDT blocks whose data is to be concatenated in the NWB file:
TDTblock_multiple = {'TremoLfpDBS-191115-100127', ...
                     'TremoLfpDBS-191115-122420', ...
                     'TremoLfpDBS-191115-125440', ...
                     'TremoLfpDBS-191115-132410', ...
                     'TremoLfpDBS-191115-140755', ...
                     'TremoLfpDBS-191115-143421'}; % manually enter in the names of the TDT blocks

% Fullpath for directory in which to save NWB file results
nwbSavePn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';


% Get time-domain related data from each TDT block to allow proper
% synchronized time data when all concatenated
nBlocks = length(TDTblock_multiple);
tic
streamStr = 'RAW8';
disp('Getting time info from all TDT blocks')
for iBlock = 1:nBlocks
    disp('Loading tdt block...')
    tdt = TDTbin2mat([tdtDataPath TDTtank '\' TDTblock_multiple{iBlock}], ...
        'STORE', streamStr, ...
        'CHANNEL', 1);
    disp('Done loading tdt block')
    %         'STORE', streamStr);
    testCh = tdt.streams.(streamStr).data(1,:); % need to get length of one chan instead of size of all chans, will flood RAM otherwise
    TDTblock_samps(iBlock) = length(testCh)
    utcStartTime{iBlock} = tdt.info.utcStartTime;
    TDTblock_dur{iBlock} = tdt.info.duration;
    fsRAW8 = tdt.streams.(streamStr).fs;

    clear tdt testCh
    
    
end
toc
disp('Done getting time info!')





% CREATE description of the RAW ElectriacalSeries object that keeps track
% of the original TDT files and how to pick out the data from each one if
% necessary. Format example:
% Concatenated TDTblocks and their sample-indices: 1) TremoLfpDBS-191115-122420[1,88621056]; ...
TDTblock_samps = TDTblock_samps';
TDTblock_sampsCumul = [ones(length(TDTblock_samps), 1), TDTblock_samps];
for iBlock = 2:nBlocks
    TDTblock_sampsCumul(iBlock,1) = TDTblock_sampsCumul(iBlock-1,2) + TDTblock_sampsCumul(iBlock,1);
    TDTblock_sampsCumul(iBlock,2) = TDTblock_sampsCumul(iBlock-1,2) + TDTblock_sampsCumul(iBlock,2);
    
end
descrString = buildElectricalSeriesDescr(TDTblock_multiple, TDTblock_sampsCumul);
           

% CREATE timestamps vector that is based on global syncronized reference
% time. If data has missing chunks of time, the time-skips will be
% reflected by discontinuities in the timestamps, according to NWB standard
% format guidelines. 

refDT = datetime(utcStartTime(1));
timeStamps = [];
for iBlock = 1:nBlocks
    % Get non-synchronized timestamps for each TDT block
    tStamps{iBlock} = (1/fsRAW8) * (0:(TDTblock_samps(iBlock)-1))';
    
    % Get synchronization reference time (seconds) for each TDT block
    refSec = seconds(datetime(utcStartTime(iBlock)) - refDT);
        
    % Re-reference4 all timestamps accordingly
    tStampsSync{iBlock} = tStamps{iBlock} + refSec;
    
    timeStamps = [timeStamps; tStampsSync{iBlock}];

end


% CREATE reference table for iteratively appending TDT data using DataPipe
% in NWB
% Fill a table with all TDTblocks (and their subdivisions) to iteratively
% append to the NWB file. Note that RAW8 is frikkin huge, so we have to use
% TDTbin2mat's T1 and T2 time limit option to chop up the data even within
% a TDT block so as not to flood my RAM...

chunkSec = 300; % divide the data into n-second chunks

TDTblock = {};
t1 = [];
t2 = [];
for iBlock = 1:nBlocks
    % Generate array of time intervals for given TDT block
    secondsDur_i = seconds(duration(TDTblock_dur{iBlock}));
    t1_i = [0:chunkSec:secondsDur_i]';
    t2_i = [t1_i(2:end); 0];
    nRows = length(t1_i);
    TDTblock_i = cell(nRows, 1);
    TDTblock_i(:) = TDTblock_multiple(iBlock);
    
    % Update the running tally of TDT block subsegments
    TDTblock = [TDTblock; TDTblock_i];
    t1 = [t1; t1_i];
    t2 = [t2; t2_i];
    
end

T_TDTload = [table(TDTblock), table(t1), table(t2)];



%% Fill out inital data in the NWB file based on the first tdt block

tdtDataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';
TDTtank = 'Uva-191115';
% TDTblock = TDTblock_multiple{1};

loadFullPn_firstBlk = [tdtDataPath TDTtank '\' TDTblock_multiple{1}];


% Initialize nwb object based on TDT metadata and user input
tdt = TDTbin2mat(loadFullPn_firstBlk, 'TYPE', {'epocs', 'snips', 'scalars'});
dt = datetime(tdt.info.date);
tt = datetime(tdt.info.utcStartTime); 
tott = dt + timeofday(tt);
sess_utcStart = datestr(tott, 'yyyy-mm-ddTHH:MM:SSZ');

% give file a globally unique id, more for machine readability than human
id = genUUID;

  

% Init NWB
nwb = NwbFile(...
    'session_start_time', sess_utcStart, ...
    'identifier', id, ...
    'session_description', 'a baseline day recording', ...
    'timestamps_reference_time', sess_utcStart);



% Fill "general" section of nwb file, constants set inside sub-function 
nwb = nwbfill_general(nwb);
nwb.general_session_id = SESSION_ID;

% Fill "electrode" information of nwb file, constants set inside sub-function
nwb = nwbfill_electrode(nwb);



%% Add RAW ephys data in to acquisition
% Store RAW8 as ElectricalSeries 

fullDataSize = [sum(TDTblock_samps), 8];


tdtDataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';
TDTtank = 'Uva-191115';



% Load Raw ephys data from tdt
streamStr = 'RAW8'; % name of the stream with raw voltage trace
tic
data = TDTbin2mat([tdtDataPath TDTtank '\' T_TDTload.TDTblock{1}], ...
    'STORE', streamStr, ...
    'T1', T_TDTload.t1(1), ...
    'T2', T_TDTload.t2(1));
% data = TDTbin2mat(loadFullPn_firstBlk, 'STORE', streamStr);
toc
rawData = data.streams.(streamStr).data'; % make m x n, where n is num channels
fsRAW8 = data.streams.RAW8.fs;
timestamps = (1/fsRAW8)*(0:(size(rawData, 1) - 1))';

% specify the electrode(s) that this signal came from according to nwb
electrodes_object_view = types.untyped.ObjectView( ...
    '/general/extracellular_ephys/electrodes');

electrode_table_region = types.hdmf_common.DynamicTableRegion( ...
    'table', electrodes_object_view, ...
    'description', 'missing', ...
    'data', (0:7)); % nwb uses zero as first index...


% Set up a DataPipe object for "data", later iteratively append data from
% more TDTblocks
rawPipe = types.untyped.DataPipe(...
    'data', rawData, ...
    'maxShape', fullDataSize, ...
    'axis', 1);


% Store data as an ElectricalSeries object
eSeriesName = 'rawVoltageTrace';
electrical_series = types.core.ElectricalSeries( ...
    'comments', 'Data stored as single-precision matlab array of 8 channels, "description" has sample indices for mapping data to original TDT blocks (if data is concatenated from multiple original files).', ...
    'data', rawPipe, ...
    'data_resolution', -1.0, ... % default for unknown
    'data_unit', 'Volts', ...
    'description', descrString, ...
    'electrodes', electrode_table_region, ... % points to value in electrode table
    'starting_time', 0.0, ... % seconds
    'starting_time_rate', fsRAW8, ...
    'timestamps', timeStamps);
    

nwb.acquisition.set(eSeriesName, electrical_series);



%% Add Startpad voltage data to Acquisition
% store StPd as TimeSeries object.
% Note: since only the first TDT block contains reach-related and
% Camera-sync related uses of the Startpad button, only the first block's
% data need be preserved in the NWB file; the rest of the TDT blocks have
% it, but nothing relevant happens


% Load Startpad data from tdt
streamStr = 'StPd'; % name of the stream with raw voltage trace
tic
data = TDTbin2mat(loadFullPn_firstBlk, 'STORE', streamStr);
toc
stPdData = data.streams.(streamStr).data'; % make m x n, where n is num channels
fsStPd = data.streams.(streamStr).fs;
% timestamps = (1/fsStPd)*(0:(size(rawData, 1) - 1))';

% % specify the electrode(s) that this signal came from according to nwb
% electrodes_object_view = types.untyped.ObjectView( ...
%     '/general/extracellular_ephys/electrodes');
% 
% electrode_table_region = types.hdmf_common.DynamicTableRegion( ...
%     'table', electrodes_object_view, ...
%     'description', 'missing', ...
%     'data', (0:7)); % nwb uses zero as first index...


% % Set up a DataPipe object for "data", later iteratively append data from
% % more TDTblocks
% rawPipe = types.untyped.DataPipe(...
%     'data', rawData, ...
%     'maxShape', fullDataSize, ...
%     'axis', 1);


% Store data as an TimeSeries object
tSeriesName = 'startPadVoltage';
time_series = types.core.TimeSeries( ...
    'comments', 'Data stored as single-precision matlab vector of voltage over time; zero V means no touch.', ...
    'data', stPdData, ...
    'data_resolution', -1.0, ... % default for unknown
    'data_unit', 'Volts', ...
    'description', ['Data from only first TDT block: ' TDTblock_multiple{1}], ...
    'starting_time', 0.0, ... % seconds
    'starting_time_rate', fsStPd);
    

nwb.acquisition.set(tSeriesName, time_series);



%% Add in IMU data stored in .csv, recroded from RaspPi

% load in csv data as table. For correcting those cases where theres 
% variables of strings-type depicting numbers (instead of just the
% numerical type for that), turns out that getting the opts object first
% automatically corrects for it...
tdt = TDTbin2mat([tdtDataPath TDTtank '\' TDTblock_multiple{1}], ...
    'TYPE', {'epocs'});

opts = detectImportOptions([tdtDataPath TDTtank '\' IMUdata]);
IMU = readtable([tdtDataPath TDTtank '\' IMUdata], opts);

% Get the detected start time of the RsP recrording in order to synchronize
% the IMU data to the NWB file's timebase (based on TDT timebase)
sync = tdt.epocs.RsP_.onset;

% Adjust the timestamp values from the IMU
IMU.TimeFromStart_s_(1:end) = IMU.TimeFromStart_s_(1:end) + sync;

% Fill in the values to be included in the NWB file
imuData = IMU{:,3:11};
time_series2 = types.core.TimeSeries( ...
    'comments', 'Data stored as double-precision columns of IMU readings over time.', ...
    'data', imuData, ...
    'data_resolution', -1.0, ... % default for unknown
    'data_unit', 'default IMU units (see IMU specs)', ...
    'description', 'AccX, AccY, AccZ, GyroX, GyroY, GyroZ, MagX, MagY, MagZ', ...
    'timestamps', IMU.TimeFromStart_s_);
    

nwb.acquisition.set('IMU', time_series2);



%% Add in timing information from the Touchscreen Task, stored within TDT
% 
% tdt = TDTbin2mat([tdtDataPath TDTtank '\' TDTblock_multiple{1}], ...
%     'TYPE', {'epocs'});






%% Add in experimental epoch information

% Original epoch interval data is in HH:MM:SS format, must convert to
% seconds for nwb:

sessionType = 'base'; % str, 'base' | 'exp'
% 
% 
% Manually enter times (only need to do this once)
% Data below can be copy-pasted as HH:MM:SS strings from <excel sheet>
session_startTime = '9:54:01';
% 
% % baseline day epochs
% naiveWashinBase = {'9:55:06', '10:56:02'}; 
% naiveTaskBase = {'10:58:06', '11:30:07'};
% naiveDbsBase = {'11:32:49', '11:48:32'};
% 
% % experiment day epochs
% naivePreExp = {};
% naiveDbsExp = {};
% harWashinExp = {};
% harTaskExp = {};
% harDbsExp = {}; 
% harWashoutExp = {};


% build a Maps object to define epochs and their time-intervals
  keySet_exp = {'naivePreExp', 'naiveDbsExp', 'harWashinExp', ...
                'harTaskExp', 'harDbsExp','harWashoutExp'};
valueSet_exp = {naivePreExp, naiveDbsExp, harWashinExp, ...
                harTaskExp, harDbsExp, harWashoutExp};

  keySet_base = {'naiveWashinBase', 'naiveTaskBase', 'naiveDbsBase'};
valueSet_base = {naiveWashinBase, naiveTaskBase, naiveDbsBase};

epcsLab = [keySet_base, keySet_exp];
vals = [valueSet_base, valueSet_exp];

% {'naiveWashinBase', 'naiveTaskBase', 'naiveDbsBase', ...
%     'naivePreExp', 'naiveDbsExp', 'harWashinExp', 'harTaskExp', ...
%     'harDbsExp','harWashoutExp'};

% vals = {naiveWashinBase, naiveTaskBase, naiveDbsBase, naivePreExp, ...
%     naiveDbsExp, harWashinExp, harTaskExp, harDbsExp, harWashoutExp};

epcs = containers.Map(epcsLab, vals);

switch sessionType
    case 'baseline'
        remove(epcs, keySet_exp);
        
    case 'exp'
        remove(epcs, keySet_base);
        
    otherwise
        error('sessionType input not recognized!')
        
end



% Now convert relevant HH:MM:SS strings into seconds (session-referenced),
% stored in "epochMap"
% epochMap = epcs;
epochKeys = keys(epcs);
nKeys = length(epochKeys);
for iKey = 1:nKeys
    epochSec{iKey} = seconds(duration(epcs(epochKeys{iKey}))) - ...
        seconds(duration(session_startTime));
    
end

epochSecMap = containers.Map(epochKeys, epochSec);
nEps = epochSecMap.Count;

% get tags
ep = epochSecMap.keys;
times = cell2mat(epochSecMap.values);
st_time = times(1:2:end);
sp_time = times(2:2:end);

% Finally, put this data into nwb Epoch Intervals section
epochs = types.core.TimeIntervals( ...
    'colnames', {'tags', 'start_time','stop_time'}, ...
    'description', 'time interval for epoch: naivePre', ...
    'id', types.hdmf_common.ElementIdentifiers( ...
        'data', 0:(nEps-1)), ...
    'start_time', types.hdmf_common.VectorData( ...
        'data', st_time, ...
        'description','start time of epoch'), ...
    'stop_time', types.hdmf_common.VectorData( ...
        'data', sp_time, ...
        'description','end time of epoch'), ...
    'tags', types.hdmf_common.VectorData( ...
        'data', ep, ...
        'description', 'Labels of separate experimental epochs, in order.'));
    
nwb.intervals_epochs = epochs;
    


%% Write the nwb file in one go

% Note: if nwb file already exists, nwbExport will error out, the origin of
% the problem being a call to hdf5lib2: "Unable to create file with 
% specified filename. Filename may have unsupported characters."
% FIRST_BLOCK = [FIRST_BLOCK 'TEST'];
tic


nwbExport(nwb, [nwbSavePn SESSION_ID '.nwb'])
% nwbExport(nwb, [NWBPATH 'TEST' '.nwb'])

disp('DONE EXPORTING NWB!');
toc

clear nwb rawData rawPipe tdt timestamps timeStamps tStamps tStampsSync


%% Read back in and test
% clear nwb
nwb2 = nwbRead([SESSION_ID '.nwb']);

nSamps = nwb2.acquisition.get(eSeriesName).data.internal.dims(1);
ch = 1;

% raw = nwb2.acquisition.get(eSeriesName).data.load([1,ch], [nSamps,ch]);
% tStamps = nwb2.acquisition.get(eSeriesName).timestamps.load;



%% Append remaining continuous data to rawVoltageTrace in acquisition group

nRows = height(T_TDTload);
streamStr = 'RAW8'; % name of the stream with raw voltage trace
disp('Appending remaining data of RAW8 to NWB file...')
tic
for iRow = 2:nRows
    disp(['Appending subsection ' num2str(iRow) ' of ' num2str(nRows)]);
    % Load Raw ephys data from tdt
    data = TDTbin2mat([tdtDataPath TDTtank '\' T_TDTload.TDTblock{iRow}], ...
    'STORE', streamStr, ...
    'T1', T_TDTload.t1(iRow), ...
    'T2', T_TDTload.t2(iRow));
    dataPart_i = data.streams.(streamStr).data'; % make m x n, where n is num channels
    nwb2.acquisition.get(eSeriesName).data.append(dataPart_i); 
    
end
disp('COMPLETE!')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
toc
% % Reset the datapipe's maxSize property from Inf to fullDataSize to
% % indicate that the dataPipe is done iterating, and ensure that no further
% % appending is possible. 
% nwb2.acquisition.get(eSeriesName).data.internal.maxSize = fullDataSize;



%% SUB-FUNCTIONS

function descrString = buildElectricalSeriesDescr(TDTblock_multiple, TDTblock_sampsCumul);
descrString = 'Concatenated TDTblocks and their sample-indices:';

nBlocks = length(TDTblock_multiple);
for iBlock = 1:nBlocks
    descrString = [descrString, ' ', num2str(iBlock), ') ', ...
        TDTblock_multiple{iBlock}, '[', ...
        num2str(TDTblock_sampsCumul(iBlock,1)), ',',...
        num2str(TDTblock_sampsCumul(iBlock,2)), '];']
    
end

end

% function nwbtable = table2DynamicTable(T)
% 
% %EXAMPLE
% %   T = table([.1, 1.5, 2.5]',[1., 2., 3.]',[0,1,0]',...
% %       'VariableNames',{'start','stop','condition'});
% %   T.Properties.Description = 'my description';
% %   T.Properties.UserData = containers.Map('source','my source');
% %nwbfile.trials = table2nwb(T)
% 
% if ismember('id', T.Properties.VariableNames)
%     id = T.id;
% else
%     id = 0:height(T)-1;
% end
% 
% nwbtable = types.core.DynamicTable( ...
%     'colnames', T.Properties.VariableNames,...
%     'description', T.Properties.Description, ...
%     'id', types.core.ElementIdentifiers('data', id));
% 
% for col = T
%     if ~strcmp(col.Properties.VariableNames{1},'id')
%         nwbtable.vectordata.set(col.Properties.VariableNames{1}, ...
%             types.core.VectorData('data', col.Variables',...
%             'description','my description'));
%     end
% end
% 
% end

function subT = getblocks(T, cols)

colNames = fieldnames(cols);
nFields = length(fieldnames(cols));

isIdx = false(height(T), nFields);
for iF = 1:nFields
    isIdx(:,iF) = strcmp(cols.(colNames{iF}), T{:,colNames{iF}});
    
end

select = all(isIdx, 2);
subT = T(select,:);


end

function [lfpConcat, tstampsConcat, fs] = concatLFP(subT, LFPLOCPATH, SUFFIX, TDTLOCPATH)


% load first block

nBlks = height(subT);
lfpConcat = [];
tstampsConcat = [];
for iBlk = 1:nBlks
    load([LFPLOCPATH subT.TDTblock{iBlk} SUFFIX])

    % Get time data using original tdt block
%     tdtdataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';
    TDTtank = [subT.TDTtank{1} '\'];
    tdtfullpath = [TDTLOCPATH TDTtank subT.TDTblock{iBlk}];
    tdt = TDTbin2mat(tdtfullpath, 'TYPE', {'scalars'});


    if iBlk == 1
        sess_tBegdt = datetime(tdt.info.utcStartTime);
    %     fs = tdt.streams.(streamField).fs; % samps/sec

    end  
    tBeg_currBlk = datetime(tdt.info.utcStartTime);


    tstamps = (1/fs) * (0:(size(lfp, 1) - 1))';

    offset = seconds(tBeg_currBlk - sess_tBegdt);
    tstampsSess = tstamps + offset;
    
    lfpConcat = [lfpConcat; lfp];
    tstampsConcat = [tstampsConcat; tstampsSess];
    
end % END forloop

end % END function

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
