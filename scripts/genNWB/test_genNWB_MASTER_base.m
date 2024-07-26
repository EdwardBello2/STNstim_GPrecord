% script for creating the "master" nwb file with all relevant metadata,
% external links to adjacent nwb files containing major data, and also
% directly write in minor data

% NOTE: if running "batch_updateAllNWB_MASTER_exp.m", must comment out the
% initial clear command and SESSION_ID hard-coded assignme in THIS script.


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

% baseline day
SESSION_ID = 'TremoLfpDBS-190923-110101';
harInjTime = '10:31:56';
harDose = '8 mg/kg';


% clear
tic
%% Get the number of samples from each of the TDT blocks to be combined

% 8 mg/kg day (1st)
% SESSION_ID = 'TremoLfpDBS-191018-100615';
% harInjTime = '10:38:43';
% harDose = '8 mg/kg';

tdtDataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';
TDTtank = ['Uva-' SESSION_ID(13:18)];
TDTblockFirst = SESSION_ID;
nwbSuffix = 'MASTER';

% Fullpath for directory in which to save NWB file results
nwbSavePn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';



%% Fill out inital data in the NWB file based on the first tdt block

% tdtDataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';
% TDTtank = 'Uva-191115';
% TDTblock = TDTblock_multiple{1};
disp('filling in general metadata...')

loadFullPn_firstBlk = [tdtDataPath TDTtank '\' TDTblockFirst];


% Initialize nwb object based on TDT metadata and user input
tdt = TDTbin2mat(loadFullPn_firstBlk, 'TYPE', {'epocs', 'snips', 'scalars'});
dt = datetime(tdt.info.date);
tt = datetime(tdt.info.utcStartTime); 
tott = dt + timeofday(tt);
sess_utcStart = datestr(tott, 'yyyy-mm-ddTHH:MM:SSZ');

% give file a globally unique id, more for machine readability than human
id = genUUID;
  
% Specify ISO 8601 extended date time string according to NWB guidelines
dateStr = util.sessionid2isodatetime(SESSION_ID);

% Init NWB
nwb = NwbFile(...
    'session_start_time', dateStr, ...
    'identifier', id, ...
    'session_description', 'a baseline day recording', ...
    'timestamps_reference_time', dateStr);



% Fill "general" section of nwb file, constants set inside sub-function 
nwb = nwbfill_general(nwb);
nwb.general_pharmacology = ['harmaline hydrochloride, dose: ' harDose];
nwb.general_notes = ['harmaline injected (IM) at time: ' harInjTime];

nwb.general_session_id = SESSION_ID;

% Fill "electrode" information of nwb file, constants set inside sub-function
nwb = nwbfill_electrode(nwb);




% %% Add Startpad voltage data to Acquisition
% % store StPd as TimeSeries object.
% % Note: since only the first TDT block contains reach-related and
% % Camera-sync related uses of the Startpad button, only the first block's
% % data need be preserved in the NWB file; the rest of the TDT blocks have
% % it, but nothing relevant happens
% 
% 
% % Load Startpad data from tdt
% streamStr = 'StPd'; % name of the stream with raw voltage trace
% tic
% data = TDTbin2mat(loadFullPn_firstBlk, 'STORE', streamStr);
% toc
% stPdData = data.streams.(streamStr).data'; % make m x n, where n is num channels
% fsStPd = data.streams.(streamStr).fs;
% 
% % Store data as an TimeSeries object
% tSeriesName = 'startPadVoltage';
% time_series = types.core.TimeSeries( ...
%     'comments', 'Data stored as single-precision matlab vector of voltage over time; zero V means no touch.', ...
%     'data', stPdData, ...
%     'data_resolution', -1.0, ... % default for unknown
%     'data_unit', 'Volts', ...
%     'description', ['Data from only first TDT block: ' TDTblockFirst], ...
%     'starting_time', 0.0, ... % seconds
%     'starting_time_rate', fsStPd);
%     
% 
% nwb.acquisition.set(tSeriesName, time_series);
% 
% disp('Done adding StartPad voltage data...')



% %% Include external link to epoch interval information
%  
% disp('linking to epoch interval data...')
% 
% intervNWBfile = [SESSION_ID '_expIntervals.nwb'];
% intervLink = types.untyped.ExternalLink([nwbSavePn intervNWBfile], ...
%     ['/intervals/epochs']);
% nwb.intervals_epochs = intervLink;
% 
% 
% 
% % %% Include external link to IMU data nwbfile
% % 
% % imuNWBfile = 'TremoLfpDBS-191115-100127_IMUdata.nwb';
% % imuLink = types.untyped.ExternalLink([nwbSavePn imuNWBfile], ...
% %     ['/acquisition/IMU']);
% % nwb.acquisition.set('IMU', imuLink);
% % disp('Done linking to IMU data...')



%% Include external link to IMU data in to acquisition

disp('linking to IMU data...')

imuNWBfile = [SESSION_ID '_IMUdata.nwb'];
ephysLink = types.untyped.ExternalLink([nwbSavePn imuNWBfile], ...
    ['/acquisition/IMU']);
nwb.acquisition.set('IMU', ephysLink);


%% Include external link to peri-reach data in to acquisition

disp('linking to periReach data...')

imuNWBfile = [SESSION_ID '_periReach.nwb'];
ephysLink = types.untyped.ExternalLink([nwbSavePn imuNWBfile], ...
    ['/acquisition/PeriReachTimes']);
nwb.acquisition.set('PeriReachTimes', ephysLink);



%% Include external link to IMU data in to processing

disp('linking to ACCproc data...')

imuNWBfile = [SESSION_ID '_ACCproc.nwb'];
ephysLink = types.untyped.ExternalLink([nwbSavePn imuNWBfile], ...
    ['/processing/ACCproc']);
nwb.processing.set('ACCproc', ephysLink);



% %% Include external link to LFP ephys data in to processing
% 
% disp('linking to processed lfp data...')
% 
% lfpNWBfile = [SESSION_ID '_proc_lfp.nwb'];
% lfpLink = types.untyped.ExternalLink([nwbSavePn lfpNWBfile], ...
%     ['/processing/ecephys']);
% nwb.processing.set('ecephys', lfpLink);
% 
% 

% %% Include external link to RAW ephys data in to acquisition
% % -------------------------------------------------------------------------
% % TEMPORARILY SUSPEND SAVING THIS, SINCE I DON'T NEED THIS DATA RIGHT
% % NOW. INCLUDING IT MAKES THE NWB FILE LOAD SUPER-SLOWLY AT THE MOMENT 
% % -------------------------------------------------------------------------
% disp('linking to raw voltage trace data...')
% 
% ephysNWBfile = [SESSION_ID '_rawEphys.nwb'];
% ephysLink = types.untyped.ExternalLink([nwbSavePn ephysNWBfile], ...
%     ['/acquisition/rawVoltageTrace']);
% nwb.acquisition.set('rawVoltageTrace', ephysLink);
% 


%% Write the nwb file in one go

% Note: if nwb file already exists, nwbExport will error out, the origin of
% the problem being a call to hdf5lib2: "Unable to create file with 
% specified filename. Filename may have unsupported characters."
% FIRST_BLOCK = [FIRST_BLOCK 'TEST'];



nwbExport(nwb, [nwbSavePn SESSION_ID '_' nwbSuffix '.nwb'])
% nwbExport(nwb, [NWBPATH 'TEST' '.nwb'])

disp('DONE EXPORTING NWB!');
toc

load handel
sound(y, Fs)



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
