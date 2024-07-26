
clear

% % 2 mg/kg day (1st)
% TDTblock = 'TremoLfpDBS-190927-100155';
% TDTtank = 'Uva-190927';
% % harInjTime = '10:31:27';
% % harDose = '2 mg/kg';

% % 4 mg/kg day (1st)
% TDTblock = 'TremoLfpDBS-191004-100637';
% TDTtank = 'Uva-191004';
% % harInjTime = '10:39:21';
% % harDose = '4 mg/kg';
 
% % 6 mg/kg day (1st)
% TDTblock = 'TremoLfpDBS-191011-104322';
% TDTtank = 'Uva-191011';
% % harInjTime = '11:14:30';
% % harDose = '6 mg/kg';
 
% % 8 mg/kg day (1st)
% TDTblock = 'TremoLfpDBS-191018-100615';
% TDTtank = 'Uva-191018';
% % harInjTime = '10:38:43';
% % harDose = '8 mg/kg';
 
% % 2 mg/kg day (2nd)
% TDTblock = 'TremoLfpDBS-191025-104651';
% TDTtank = 'Uva-191025';
% % harInjTime = '11:17:58';
% % harDose = '2 mg/kg';
 
% 4 mg/kg day (2nd)
TDTblock = 'TremoLfpDBS-191101-101430';
TDTtank = 'Uva-191101';
% harInjTime = '10:48:15';
% harDose = '4 mg/kg';
 
% % 6 mg/kg day (2nd)
% TDTblock = 'TremoLfpDBS-191108-101829';
% TDTtank = 'Uva-191108';
% % harInjTime = '10:50:09';
% % harDose = '6 mg/kg';
 
% % 8 mg/kg day (2nd)
% TDTblock = 'TremoLfpDBS-191115-100127';
% TDTtank = 'Uva-191115';
% % harInjTime = '10:31:56';
% % harDose = '8 mg/kg';


% baseline day
TDTblock = 'TremoLfpDBS-190923-110101';
TDTtank = 'Uva-190923';


%% 


% MANUAL INPUTS:
tdtDataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';
% TDTtank = 'Uva-190927';
% IMUdata = 'TremoLfpDBS-191115-100127_IMU.csv';
% TDTblock = 'TremoLfpDBS-190927-100155';
SESSION_ID = TDTblock;

% Fullpath for directory in which to save NWB file results
nwbSavePn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';
nwbSuffix = 'periReach';


%% Calculate Reach-related timings

% Load tdt file
DataLocalPn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';
% TDTtank = 'Uva-191025';
% TDTblock = 'TremoLfpDBS-191025-104651';

tdt = TDTbin2mat([DataLocalPn TDTtank '\' TDTblock], ...
    'TYPE', {'epocs'});

tdt2 = TDTbin2mat([DataLocalPn TDTtank '\' TDTblock], ...
    'STORE', 'StPd');

stpdData = tdt2.streams.StPd.data;
fs = tdt2.streams.StPd.fs;
stpdTime = (1/fs) * (0:(length(stpdData)-1)) + tdt2.streams.StPd.startTime;
stpdMidpt = (max(stpdData) + min(stpdData)) / 2;



% Gather data for SUCCESSFUL REACHES only

% Get ReachCompletion times
reachEnd = tdt.epocs.PC2_.onset;


% Get Target presentation onset times
targSpanAll = [tdt.epocs.PC1_.onset tdt.epocs.PC1_.offset];
nReaches = length(reachEnd);
targSpanReach = zeros(nReaches, 2);

% Get the target presentation span for each successful reach completion
for iReach = 1:nReaches
    % check each reach completion time for the immediately preceding target
    % span
    iRchTime = reachEnd(iReach) - 0.05;
    
    isInReach = (targSpanAll(:,1) < iRchTime) & (targSpanAll(:,2) > iRchTime);
    if sum(isInReach) ~= 1, error('aaaaa'); end
    
    targSpanReach(iReach,1:2) = targSpanAll(isInReach,1:2);
       
end
targOnset = targSpanReach(:,1);


% Get Reach Onset time (i.e. time when paw leaves button; stpd voltage
% drops)
reachOnset = zeros(nReaches, 1);
for iReach = 1:nReaches
    % get window of stpd data within reach target span
    isInSpan = (stpdTime >= targSpanReach(iReach,1)) & ...
        (stpdTime < targSpanReach(iReach,2));
    stpdTimeWin = stpdTime(isInSpan);
    stpdDataWin = stpdData(isInSpan);
    
    % Find first timepoint where stpd voltage falls below the threshold
    idxBelow = find((stpdDataWin < stpdMidpt), 1);
    reachOnset(iReach,1) = stpdTimeWin(idxBelow);
    
end



%% Create NWB object

% loadFullPn_firstBlk = [tdtDataPath TDTtank '\' TDTblock];
% 
% % Initialize nwb object based on TDT metadata and user input
% tdt = TDTbin2mat(loadFullPn_firstBlk, 'TYPE', {'epocs', 'snips', 'scalars'});
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



% %% Add in experimental epoch information
% 
% % Original epoch interval data is in HH:MM:SS format, must convert to
% % seconds for nwb:
% 
% sessionType = 'exp'; % str, 'base' | 'exp'
% % 
% % 
% % Manually enter times (only need to do this once)
% % Data below can be copy-pasted as HH:MM:SS strings from <excel sheet>
% session_startTime = '10:01:58';
% 
% % baseline day epochs; may just be dummy values if 'exp' the session type
% naiveWashinBase = {'::', '::'}; 
% naiveTaskBase = {'::', '::'};
% naiveDbsBase = {'::', '::'};
% 
% % experiment day epochs
%   naivePreExp = {'10:04:07', '10:12:36'};
%   naiveDbsExp = {'10:13:09', '10:28:41'};
%  harWashinExp = {'10:31:33', '11:31:43'};
%    harTaskExp = {'11:33:27', '12:04:17'};
%     harDbsExp = {'12:06:20', '12:22:06'}; 
% harWashoutExp = {'12:23:42', '14:39:34'};
% 
% 
% % build a Maps object to define epochs and their time-intervals
%   keySet_exp = {'naivePreExp', 'naiveDbsExp', 'harWashinExp', ...
%                 'harTaskExp', 'harDbsExp','harWashoutExp'};
% valueSet_exp = {naivePreExp, naiveDbsExp, harWashinExp, ...
%                 harTaskExp, harDbsExp, harWashoutExp};
% 
%   keySet_base = {'naiveWashinBase', 'naiveTaskBase', 'naiveDbsBase'};
% valueSet_base = {naiveWashinBase, naiveTaskBase, naiveDbsBase};
% 
% epcsLab = [keySet_base, keySet_exp];
% vals = [valueSet_base, valueSet_exp];
% 
% % {'naiveWashinBase', 'naiveTaskBase', 'naiveDbsBase', ...
% %     'naivePreExp', 'naiveDbsExp', 'harWashinExp', 'harTaskExp', ...
% %     'harDbsExp','harWashoutExp'};
% 
% % vals = {naiveWashinBase, naiveTaskBase, naiveDbsBase, naivePreExp, ...
% %     naiveDbsExp, harWashinExp, harTaskExp, harDbsExp, harWashoutExp};
% 
% epcs = containers.Map(epcsLab, vals);
% 
% switch sessionType
%     case 'baseline'
%         remove(epcs, keySet_exp);
%         
%     case 'exp'
%         remove(epcs, keySet_base);
%         
%     otherwise
%         error('sessionType input not recognized!')
%         
% end
% 
% 
% 
% % Now convert relevant HH:MM:SS strings into seconds (session-referenced),
% % stored in "epochMap"
% % epochMap = epcs;
% epochKeys = keys(epcs);
% nKeys = length(epochKeys);
% for iKey = 1:nKeys
%     epochSec{iKey} = seconds(duration(epcs(epochKeys{iKey}))) - ...
%         seconds(duration(session_startTime));
%     
% end
% 
% epochSecMap = containers.Map(epochKeys, epochSec);
% nEps = epochSecMap.Count;
% 
% % get tags
% ep = epochSecMap.keys;
% times = cell2mat(epochSecMap.values);
% st_time = times(1:2:end);
% sp_time = times(2:2:end);
% 
% tagidx = [4 3 2 5 1 0]; % indices to indicate the order of occurrence of these intervals, they don't get put in order
% 
% 
% 
% %% Fill NWB Behavioral Events object with the data
% 
% % Finally, put this data into nwb Epoch Intervals section
% epochs = types.core.TimeIntervals( ...
%     'colnames', {'tags', 'start_time','stop_time'}, ...
%     'description', 'time intervals for experimental epochs', ...
%     'id', types.hdmf_common.ElementIdentifiers( ...
%         'data', tagidx), ...
%     'start_time', types.hdmf_common.VectorData( ...
%         'data', st_time, ...
%         'description','start time of epoch'), ...
%     'stop_time', types.hdmf_common.VectorData( ...
%         'data', sp_time, ...
%         'description','end time of epoch'), ...
%     'tags', types.hdmf_common.VectorData( ...
%         'data', ep, ...
%         'description', 'Labels of separate experimental epochs, may not be in chron order.'));
%     
% nwb.intervals_epochs = epochs;
%     
% disp('Done adding experimental intervals...')
% 
% 
% 
% % Fill in the values to be included in the NWB file
% reachData = [targOnset, reachOnset, reachEnd];
% time_series2 = types.core.TimeSeries( ...
%     'comments', 'Data stored as double-precision columns of Reach-related event times.', ...
%     'data', reachData, ...
%     'data_resolution', -1.0, ... % default for unknown
%     'data_unit', 'seconds', ...
%     'description', 'TargetOnset, ReachOnset, ReachEnd', ...
%     'timestamps', IMU.TimeFromStart_s_, ...
%     'control', ctrl, ...
%     'control_description', 'numbers correspond to sample-by-sample acquisition attempts from recording system');
%     

% nwb.acquisition.set('ReachData', time_series2);

% target onset
data = 1:length(targOnset);
targetOnsetTseries = types.core.TimeSeries( ...
    'comments', 'Number is the number of reaches.', ...
    'data', data, ...
    'data_unit', 'reach-count', ...
    'data_resolution', -1.0, ... % default for unknown
    'description', 'TargetOnset, ReachOnset, ReachEnd', ...
    'timestamps', targOnset);
    
% reach onset
data = 1:length(reachOnset);
reachOnsetTseries = types.core.TimeSeries( ...
    'comments', 'Number is the number of reaches.', ...
    'data', data, ...
    'data_unit', 'reach-count', ...
    'data_resolution', -1.0, ... % default for unknown
    'description', 'TargetOnset, ReachOnset, ReachEnd', ...
    'timestamps', reachOnset);

% reach end time
data = 1:length(reachEnd);
reachEndTseries = types.core.TimeSeries( ...
    'comments', 'Number is the number of reaches.', ...
    'data', data, ...
    'data_unit', 'reach-count', ...
    'data_resolution', -1.0, ... % default for unknown
    'description', 'TargetOnset, ReachOnset, ReachEnd', ...
    'timestamps', reachEnd);


BehavioralEvents = types.core.BehavioralEvents;
BehavioralEvents.timeseries.set('targetOnset', targetOnsetTseries);
BehavioralEvents.timeseries.set('reachOnset', reachOnsetTseries);
BehavioralEvents.timeseries.set('reachEnd', reachEndTseries);

nwb.acquisition.set('PeriReachTimes', BehavioralEvents);
%% Write the nwb file in one go

% Note: if nwb file already exists, nwbExport will error out, the origin of
% the problem being a call to hdf5lib2: "Unable to create file with 
% specified filename. Filename may have unsupported characters."
% FIRST_BLOCK = [FIRST_BLOCK 'TEST'];
tic


nwbExport(nwb, [nwbSavePn SESSION_ID '_' nwbSuffix '.nwb'])
% nwbExport(nwb, [NWBPATH 'TEST' '.nwb'])

disp('DONE EXPORTING NWB!');
toc



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
