% % 2 mg/kg day (1st)
% SESSION_ID = 'TremoLfpDBS-190927-100155';
% harInjTime = '10:31:27';
% harDose = '2 mg/kg';

% % 4 mg/kg day (1st)
% SESSION_ID = 'TremoLfpDBS-191004-100637';
% harInjTime = '10:39:21';
% harDose = '4 mg/kg';

% % 6 mg/kg day (1st)
% SESSION_ID = 'TremoLfpDBS-191011-104322';
% harInjTime = '11:14:30';
% harDose = '6 mg/kg';

% % 8 mg/kg day (1st)
% SESSION_ID = 'TremoLfpDBS-191018-100615';
% harInjTime = '10:38:43';
% harDose = '8 mg/kg';

% % 2 mg/kg day (2nd)
% SESSION_ID = 'TremoLfpDBS-191025-104651';
% harInjTime = '11:17:58';
% harDose = '2 mg/kg';

% % 4 mg/kg day (2nd)
% SESSION_ID = 'TremoLfpDBS-191101-101430';
% harInjTime = '10:48:15';
% harDose = '4 mg/kg';

% 6 mg/kg day (2nd)
SESSION_ID = 'TremoLfpDBS-191108-101829';
harInjTime = '10:50:09';
harDose = '6 mg/kg';

% % 8 mg/kg day (2nd)
% SESSION_ID = 'TremoLfpDBS-191115-100127';
% harInjTime = '10:31:56';
% harDose = '8 mg/kg';



%% 
% load in csv data as table. For correcting those cases where theres 
% variables of strings-type depicting numbers (instead of just the
% numerical type for that), turns out that getting the opts object first
% automatically corrects for it...
% clear

% MANUAL INPUTS:
tdtDataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';
projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';
IMUdataPath = [projRootPath 'Data Acquisition\Uva dose study\IMUdata\'];
% SESSION_ID = 'TremoLfpDBS-190927-100155';

IMUdata = [SESSION_ID '_IMU.mat'];
TDTblock = SESSION_ID;
TDTtank = ['Uva-' SESSION_ID(13:18)];



% Fullpath for directory in which to save NWB file results
nwbSavePn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';
nwbSuffix = 'IMUdata';


%% Create NWB object

loadFullPn_firstBlk = [tdtDataPath TDTtank '\' TDTblock];

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



%%   Add in IMU data recorded from RaspPi

% data is stored in matlab tables, after having been processed to correct
% any errors during acquisition

% tdt = TDTbin2mat([tdtDataPath TDTtank '\' TDTblock], ...
%     'TYPE', {'epocs'});

% load matlab table(s) with IMU data and timestamps, also sync info
load([IMUdataPath IMUdata], 'IMUcell');

nTabs = numel(IMUcell);
IMU = []; % initialize final table
ctrl = []; % tracking the sample numbers for each segment of recordings that day
for iTab = 1:nTabs
    tempIMU = IMUcell{iTab};
    % Get the detected start time of the RsP recrording in order to synchronize
    % the IMU data to the NWB file's timebase (based on TDT timebase)
    syncTDT = tempIMU.Properties.UserData.sync;

    % Adjust the timestamp values from the IMU
    tempIMU.TimeFromStart_s_(1:end) = tempIMU.TimeFromStart_s_(1:end) + syncTDT;
    i_acq = [0:(height(tempIMU)-1)]';

    
    IMU = [IMU; tempIMU];
    ctrl = [ctrl; (iTab * ones(length(i_acq), 1))];
    
    
end


% Fill in the values to be included in the NWB file
imuData = IMU{:,2:10};
time_series2 = types.core.TimeSeries( ...
    'comments', 'Data stored as double-precision columns of IMU readings over time.', ...
    'data', imuData, ...
    'data_resolution', -1.0, ... % default for unknown
    'data_unit', 'default IMU units (see IMU specs)', ...
    'description', 'AccX, AccY, AccZ, GyroX, GyroY, GyroZ, MagX, MagY, MagZ', ...
    'timestamps', IMU.TimeFromStart_s_, ...
    'control', ctrl, ...
    'control_description', 'numbers correspond to sample-by-sample acquisition attempts from recording system');
    

nwb.acquisition.set('IMU', time_series2);



%% Write the nwb file in one go

% Note: if nwb file already exists, nwbExport will error out, the origin of
% the problem being a call to hdf5lib2: "Unable to create file with 
% specified filename. Filename may have unsupported characters."
% FIRST_BLOCK = [FIRST_BLOCK 'TEST'];
tic

fullSavePath = [nwbSavePn SESSION_ID '_' nwbSuffix '.nwb'];

% Next two lines essentially let you overwrite an existing version of this
% nwb file
nwbnrtl.io.checkdelete(fullSavePath); 
nwbExport(nwb, fullSavePath)
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
