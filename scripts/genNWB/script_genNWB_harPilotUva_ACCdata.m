% script for converting original Mela accel data into NWB format. Done so that
% data between Mela, Kiwi, and Uva may be easily compared and the same
% scripts run on them all, since original formats are all frikkiin
% different. 

%% NWB code for converting original accel data to NWB format

% Collected on microstrain G-link hardware
clear; close all

% % G-link sampling freq when streaming 3 channels (hardware design)
% fs = 617; % samples/second


% MANUAL INPUTS:
projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';
dataAcqPn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\';
tdtPn = [dataAcqPn 'TDTdata\'];
% ACCdataPn = [dataAcqPn 'Uva\Harmaline\'];
SESSION_ID_prefix = 'harPilotUva-';


% Fullpath for directory in which to save NWB file results
nwbSavePn = [dataAcqPn 'NWBdata\UvaPilot\'];
nwbSuffix = '_ACCdata';


% Speicfy which day's session-data is being concatenated and stored in NWB
% format
sessionNum = 3;


% INPUT METADATA TABLES
% Read in metadata table for Uva's recordings
tablePath = 'Data Acquisition\Uva Pilot\';
metaTab = readtable([projRootPath tablePath 'acquisitionMetadata_Uva.csv']);
harTab = readtable([projRootPath tablePath 'harmalineSessionData.csv']);
fullMeta = join(metaTab, harTab);

% add in data for time since harmaline injection
harRefTime = fullMeta.startTime - fullMeta.harInjTime;
fullMeta = [fullMeta, table(harRefTime)];

% Remove rows not pertaining to the experiment
isBad = strcmp(fullMeta.recComment, 'not part of experiment') | ...
    strcmp(fullMeta.recComment, 'fake human tremor');
fullMeta = fullMeta(~isBad,:);



%%   Collect Accel data recorded from GLink system for the entire session

SessN = fullMeta(fullMeta.session == sessionNum,:);
infmt = 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z';
dstr = SessN.dateStr{1};
refDatetime = datetime(dstr, 'InputFormat', infmt);

% get full session id label for naming NWB file based on first file
SESSION_ID = [SESSION_ID_prefix,  dstr(3:4), dstr(6:7), dstr(9:10), ...
    '-', dstr(12:13) dstr(15:16) dstr(18:19)]

% CREATE timestamps vector that is based on global syncronized reference
% time. If data has missing chunks of time, the time-skips will be
% reflected by discontinuities in the timestamps, according to NWB standard
% format guidelines. 
timeStamps = [];
acc = [];
acqAttempt = [];
ctrl = [];
nFiles = height(SessN);
for iFile = 1:nFiles
    
    % Get non-synchronized timestamps for each file in session
    tdt = TDTbin2mat([tdtPn SessN.TDTtank{iFile} '\' SessN.TDTblock{iFile}], ...
        'STORE', 'Acc3');


    [i_tst, i_acc, i_acq] = getAccInfoUva(tdt);

    % Get synchronization reference time (seconds) for each TDT block
    refSec = seconds(datetime(SessN.dateStr{iFile}, 'InputFormat', infmt) - refDatetime);
        
    % Re-reference4 all timestamps accordingly
    i_tstSync = i_tst + refSec;
    
    timeStamps = [timeStamps; i_tstSync];
    acc = [acc; i_acc];
    acqAttempt = [acqAttempt; i_acq];
    
    ctrl = [ctrl; (iFile * ones(length(i_tst), 1))];
    
end



%% Create NWB object and Add in data to NWB object and save final result
% This file to be accessed by the "external link" property of the MASTER
% nwb file

% give file a globally unique id, more for machine readability than human
id = genUUID;

% Init NWB
nwb = NwbFile(...
    'session_start_time', SessN.dateStr{1}, ...
    'identifier', id, ...
    'session_description', 'a harmaline-recording pilot study', ...
    'timestamps_reference_time', SessN.dateStr{1});

% Fill in the values to be included in the NWB file
% imuData = IMU{:,2:10};
time_series2 = types.core.TimeSeries( ...
    'comments', 'Data stored as double-precision columns of IMU readings over time.', ...
    'data', acc, ...
    'data_resolution', -1.0, ... % default for unknown
    'data_unit', 'missing', ...
    'description', 'AccX, AccY, AccZ', ...
    'timestamps', timeStamps, ...
    'starting_time', 0, ...
    'starting_time_rate', 1.0173e+03, ...
    'control', ctrl, ...
    'control_description', 'numbers label the data samples correspond to N different recordings concatenated within this NWB file');
    

nwb.acquisition.set('ACC', time_series2);



% Write the nwb file in one go
tic

fullSavePath = [nwbSavePn SESSION_ID nwbSuffix '.nwb'];

% Next two lines essentially let you overwrite an existing version of this
% nwb file
nwbnrtl.io.checkdelete(fullSavePath); 
nwbExport(nwb, fullSavePath)

% nwbExport(nwb, [NWBPATH 'TEST' '.nwb'])

disp('DONE EXPORTING NWB!');
toc



%% Test loading and navigating the new nwb file

nwb2 = nwbRead([nwbSavePn SESSION_ID nwbSuffix '.nwb'])



%% SUB-FUNCTIONS

function [i_tst, i_acc, i_acq] = getAccInfoUva(tdt)

fs = tdt.streams.Acc3.fs; % hz
i_acc = double(tdt.streams.Acc3.data');
sampNum = length(i_acc);
i_acq = [0:(sampNum-1)]';
i_tst = (1/fs) * i_acq;
    
end

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


