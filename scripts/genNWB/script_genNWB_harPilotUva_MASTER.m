% script for creating the "master" nwb file with all relevant metadata,
% external links to adjacent nwb files containing major data, and also
% directly write in minor data

% NOTE: if running "batch_updateAllNWB_MASTER_exp.m", must comment out the
% initial clear command and SESSION_ID hard-coded assignme in THIS script.

% Day 1
% SESSION_ID = 'harPilotUva-180219-132618';
% harInjTime = '13:45:00';
% harDose = '8 mg/kg';

% Day 2
% SESSION_ID = 'harPilotUva-180221-113120';
% harInjTime = '11:30:00';
% harDose = '12 mg/kg';

% Day 3
% SESSION_ID = 'harPilotUva-180305-103016';
% harInjTime = '10:53:00';
% harDose = '12 mg/kg';

tic
%% Get the number of samples from each of the TDT blocks to be combined
clear

SESSION_ID = 'harPilotUva-180221-113120';
harInjTime = '11:30:00';
harDose = '12 mg/kg';

tdtDataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';
TDTtank = ['Uva-' SESSION_ID(13:18)];
TDTblockFirst = SESSION_ID;
nwbSuffix = '_MASTER';

% Fullpath for directory in which to save NWB file results
nwbSavePn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';
nwbFullSavePn = [nwbSavePn 'UvaPilot\'];



%% Fill out inital data in the NWB file based on the first tdt block

disp('filling in general metadata...')


% Specify ISO 8601 extended date time string according to NWB guidelines
dateStr = util.sessionid2isodatetime(SESSION_ID);

% give file a globally unique id, more for machine readability than human
id = genUUID;

  

% Init NWB
nwb = NwbFile(...
    'session_start_time', dateStr, ...
    'identifier', id, ...
    'session_description', 'a harmaline pilot day recording', ...
    'timestamps_reference_time', dateStr);



% Fill "general" section of nwb file, constants set inside sub-function 
% nwb = nwbfill_general(nwb);
nwb.general_session_id = SESSION_ID;
nwb.general_pharmacology = ['harmaline hydrochloride, dose: ' harDose];
nwb.general_notes = ['harmaline injected (IM) at time: ' harInjTime];

% Fill "electrode" information of nwb file, constants set inside sub-function
% nwb = nwbfill_electrode(nwb);



%% Include external link to IMU data in to acquisition

disp('linking to ACC data...')

imuNWBfile = [SESSION_ID '_ACCdata.nwb'];
ephysLink = types.untyped.ExternalLink([nwbFullSavePn imuNWBfile], ...
    ['/acquisition/ACC']);
nwb.acquisition.set('ACC', ephysLink);



%% Include external link to IMU data in to processing

disp('linking to ACCproc data...')

imuNWBfile = [SESSION_ID '_ACCproc.nwb'];
ephysLink = types.untyped.ExternalLink([nwbFullSavePn imuNWBfile], ...
    ['/processing/ACCproc']);
nwb.processing.set('ACCproc', ephysLink);



%% Include external link to movement artifact interva data in to acquisition

disp('linking to movement artifact interval data...')

movNWBfile = [SESSION_ID '_moveArt.nwb'];
ephysLink = types.untyped.ExternalLink([nwbFullSavePn movNWBfile], ...
    ['/acquisition/movementArtifacts']);
nwb.acquisition.set('movementArtifacts', ephysLink);



%% Write the nwb file in one go

% Note: if nwb file already exists, nwbExport will error out, the origin of
% the problem being a call to hdf5lib2: "Unable to create file with 
% specified filename. Filename may have unsupported characters."
% FIRST_BLOCK = [FIRST_BLOCK 'TEST'];



nwbExport(nwb, [nwbFullSavePn SESSION_ID nwbSuffix '.nwb'])
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
