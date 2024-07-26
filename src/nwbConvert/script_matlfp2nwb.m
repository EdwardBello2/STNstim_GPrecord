% Convert all relevant data from TDT blocks from entire recording session
% into one nwb file
%
%

% Authro: Ed Bello
% Created: 4/7/2020
%
%
% ChangeLog:
% 
%
% To-Do
% - use UUID for globally unique file identifier
% - implement record for which TDT files went into this NWB file
%
%
%% Specify paths needed for inputs and outputs

% Add (temporarily) folder & subfolders of Git version-controlled code for Harmaline project
addpath(genpath('C:\Users\bello043\Documents\GitHub\Thalamic-DBS-HarmalineTremor'));

% Local path to google drive project folder
PROJROOTPATH = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';

% Local path to Original TDT data
TDTLOCPATH = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\'; 

% Local path to LFP data
LFPLOCPATH = [PROJROOTPATH 'Data Processing\LFPdownsamp\'];
SUFFIX = '_diffLFP'; % append to TDTblock filename to get .mat file containing converted LFP data

% Local NWB output 
NWBPATH = [PROJROOTPATH 'Data Processing\NWB\'];


%% CONSTANTS & PARAMETERS

NHP_RARID = '13LP1';

% Read in table with metadata on all TDT blocks (harmaline details, etc.)
T = readtable([LFPLOCPATH 'TDTblockExpMetadata.xlsx']);

% get subselection from table for a given session
cols.harmalineDose = '8 mg/kg'; % str, '2 mg/kg' | '4 mg/kg' | '6 mg/kg' | '8 mg/kg'
cols.exposure = '2nd'; % str, '1st' | '2nd'
cols.sessionType = 'exp'; % str, 'baseline' | 'exp'
subT = getblocks(T, cols);


% % For Epoch time intervals:
% sessionType = 'base'; % str, 'base' | 'exp'

% Manually enter times (only need to do this once)
% Data below can be copy-pasted as HH:MM:SS strings from <excel sheet>
session_startTime = '10:01:31';
session_stopTime = '15:03:39';

% session_startTime = '::';
% session_stopTime = '::';


% baseline day epochs
naiveWashinBase = {'::', '::'}; % {'::', '::'}
naiveTaskBase = {'::', '::'}; % {'::', '::'}
naiveDbsBase = {'::', '::'}; % {'::', '::'}


% experiment day epochs
naivePreExp = {session_startTime, '10:12:34'};
naiveDbsExp = {'10:12:34', '10:28:39'};
harWashinExp = {'10:33:00', '11:33:04'};
harTaskExp = {'11:34:33', '12:05:01'};
harDbsExp = {'12:06:28', '12:23:40'};
harWashoutExp = {'12:24:23', session_stopTime};

% naivePreExp = {session_startTime, '::'};
% naiveDbsExp = {'::', '::'};
% harWashinExp = {'::', '::'};
% harTaskExp = {'::', '::'};
% harDbsExp = {'::', '::'};
% harWashoutExp = {'::', session_stopTime};




% Concatenate LFP data and session-referenced timestamps from subT
[lfpData, timestamps, fs] = concatLFP(subT, LFPLOCPATH, SUFFIX, TDTLOCPATH);



%% Initialize nwb object based on TDT metadata and user input

FIRST_BLOCK = subT.TDTblock{1};
TDTtank = [subT.TDTtank{1} '\'];
fulltdtpath = [TDTLOCPATH TDTtank FIRST_BLOCK];
tdt = TDTbin2mat(fulltdtpath, 'TYPE', {'epocs', 'snips', 'scalars'});
dt = datetime(tdt.info.date);
tt = datetime(tdt.info.utcStartTime); 
tott = dt + timeofday(tt);
sess_utcStart = datestr(tott, 'yyyy-mm-ddTHH:MM:SSZ');
% id = [sess_utcStart(1:4), sess_utcStart(6:7), sess_utcStart(9:10), '_', ...
%       NHP_RARID];

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
nwb.general_session_id = FIRST_BLOCK;

% Fill "electrode" information of nwb file, constants set inside sub-function
nwb = nwbfill_electrode(nwb);



%% Add LFP data in to processing
% Store LFP (generally downsampled and/or filtered data) as an 
% ElectricalSeries in a processing module called 'ecephys'.


% specify the electrode(s) that this signal came from according to nwb
electrodes_object_view = types.untyped.ObjectView( ...
    '/general/extracellular_ephys/electrodes');

electrode_table_region = types.hdmf_common.DynamicTableRegion( ...
    'table', electrodes_object_view, ...
    'description', 'missing', ...
    'data', (0:7)); % nwb uses zero as first index...


% store lfps as electrical series object in nwb under Processing
electrical_series = types.core.ElectricalSeries( ...
    'comments', 'data stored as single-precision matlab array, indices for mapping current data to previous TDT blocks is in "description"', ...
    'data', lfpData, ...
    'data_resolution', -1.0, ... % default for unknown
    'data_unit', 'Volts', ...
    'description', 'missing', ...
    'electrodes', electrode_table_region, ... % points to value in electrode table
    'starting_time', 0.0, ... % seconds
    'starting_time_rate', fs, ...
    'timestamps', timestamps);
    
    

ecephys_module = types.core.ProcessingModule(...
    'description', 'holds extracellular electrophysiology data');
ecephys_module.nwbdatainterface.set('LFP', ...
    types.core.LFP('lfp', electrical_series));
nwb.processing.set('ecephys', ecephys_module);



%% Add in experimental epoch information

% Original epoch interval data is in HH:MM:SS format, must convert to
% seconds for nwb:

% sessionType = 'base'; % str, 'base' | 'exp'
% 
% 
% % Manually enter times (only need to do this once)
% % Data below can be copy-pasted as HH:MM:SS strings from <excel sheet>
% session_startTime = '9:54:01';
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

switch cols.sessionType
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
    
    
% 
% 
% 
% 
% % borrowed code:
% % Re-reference clock times of epochs to session_strattime
%     Tref = subT;
%     refCol = 5; % column for sessionStartType duration object
%     Tref{:,5:end} = subT{:,5:end} - subT{:,refCol};
% 
%     harDeliv = seconds(Tref.harDelivered(:));
%     harDeliv = harDeliv/60; % convert to minutes
% 
%     % pre-harmaline naive epoch
%     t_naivePre = [seconds(Tref.session_startTime(:)), seconds(Tref.naivePre_end(:))];
% 
%     % washin epoch
%     t_washin = [seconds(Tref.harDelivered(:)), seconds(Tref.washin_end(:))];
% 
%     % washout epoch
%     t_washout = [seconds(Tref.washout_beg(:)), seconds(Tref.session_endTime(:))];
% 
% 
%     % Get indices in data
%     isNaivePre = (tst >= t_naivePre(1)) & (tst < t_naivePre(2));
%       isWashin = (tst >= t_washin(1))   & (tst < t_washin(2));
%      isWashout = (tst >= t_washout(1))  & (tst < t_washout(2));
% 
%     isWashExp = isNaivePre | isWashin | isWashout;
% 
%     
%     % Get sub-selection of LFP data and timestamps
%     lfp_naivePre = lfpFilt(isNaivePre,cpair);
%     tst_naivePre = tst(isNaivePre);
% 
%     lfp_washin = lfpFilt(isWashin,cpair);
%     tst_washin = tst(isWashin);
% 
%     lfp_washout = lfpFilt(isWashout,cpair);
%     tst_washout = tst(isWashout);
%     
    
    
    

%% Write the nwb file in one go

% Note: if nwb file already exists, nwbExport will error out, the origin of
% the problem being a call to hdf5lib2: "Unable to create file with 
% specified filename. Filename may have unsupported characters."
% FIRST_BLOCK = [FIRST_BLOCK 'TEST'];
nwbExport(nwb, [NWBPATH FIRST_BLOCK '.nwb'])
% nwbExport(nwb, [NWBPATH 'TEST' '.nwb'])

disp('DONE EXPORTING NWB!');



%% SUB-FUNCTIONS

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
