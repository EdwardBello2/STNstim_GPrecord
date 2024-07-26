% script that makes all 8 LFP process NWB files. Based on simply reapeating
% the code "script_ --- " of same name, for 8 different SESSION ID's,
% basicially automate myself so taht the 8 files get made in 40 mins

% 
% clear
% tic
% % Collect sample lims for the multiple tdt "files" represented in NWB
% % all values to be found in the string within raw.description
% sampLims = [1, 414740480; ...
%             414740481, 502513664; ...
%             502513665, 590540800; ...
%             590540801, 676052992; ...
%             676052993, 764444672; ...
%             764444673, 811200512];
% func_genNWB_procLFP('TremoLfpDBS-190927-100155', sampLims)
% toc
% disp('done with 1 of 8')


clear
tic
% Collect sample lims for the multiple tdt "files" represented in NWB
% all values to be found in the string within raw.description
sampLims = [1, 439635968; ...
            439635969, 528179200; ...
            528179201, 618979328; ...
            618979329, 710524928; ...
            710524929, 798605312; ...
            798605313, 908722176];
func_genNWB_procLFP('TremoLfpDBS-191004-100637', sampLims)
toc
disp('done with 2 of 8')


clear
tic
% Collect sample lims for the multiple tdt "files" represented in NWB
% all values to be found in the string within raw.description
sampLims = [1, 434978816; ...
            434978817, 524599296; ...
            524599297, 616865792; ...
            616865793, 714428416; ...
            714428417, 805736448; ...
            805736449, 918175744];
func_genNWB_procLFP('TremoLfpDBS-191011-104322', sampLims)
toc
disp('done with 3 of 8')


clear
tic
% Collect sample lims for the multiple tdt "files" represented in NWB
% all values to be found in the string within raw.description
sampLims = [1, 423178240; ...
            423178241, 511459328; ...
            511459329, 601247744; ...
            601247745, 691814400; ...
            691814401, 763494400; ...
            763494401, 853143552];
func_genNWB_procLFP('TremoLfpDBS-191018-100615', sampLims)
toc
disp('done with 4 of 8')


clear
tic
% Collect sample lims for the multiple tdt "files" represented in NWB
% all values to be found in the string within raw.description
sampLims = [1, 422174720; ...
            422174721, 517246976; ...
            517246977, 605245440; ...
            605245441, 694071296; ...
            694071297, 782516224; ...
            782516225, 870125568];
func_genNWB_procLFP('TremoLfpDBS-191025-104651', sampLims)
toc
disp('done with 5 of 8')


clear
tic
% Collect sample lims for the multiple tdt "files" represented in NWB
% all values to be found in the string within raw.description
sampLims = [1, 427589632; ...
            427589633, 515145728; ...
            515145729, 623554560; ...
            623554561, 701145088; ...
            701145089, 789983232; ...
            789983233, 882806784];
func_genNWB_procLFP('TremoLfpDBS-191101-101430', sampLims)
toc
disp('done with 6 of 8')


clear
tic
% Collect sample lims for the multiple tdt "files" represented in NWB
% all values to be found in the string within raw.description
sampLims = [1, 427982848; ...
            427982849, 512856064; ...
            512856065, 603758592; ...
            603758593, 692887552; ...
            692887553, 782954496; ...
            782954497, 870563840];
func_genNWB_procLFP('TremoLfpDBS-191108-101829', sampLims)
toc
disp('done with 7 of 8')


clear
tic
% Collect sample lims for the multiple tdt "files" represented in NWB
% all values to be found in the string within raw.description
sampLims = [1, 418291712; ...
            418291713, 506912768; ...
            506912769, 592953344; ...
            592953345, 720818176; ...
            720818177, 798003200; ...
            798003201, 883650560];
func_genNWB_procLFP('TremoLfpDBS-191115-100127', sampLims)
toc
disp('done with 8 of 8')



function func_genNWB_procLFP(SESSION_ID, sampLims)

%% CONSTANTS

% SESSION_ID = 'TremoLfpDBS-190927-100155';

% % bipolar electrode pair:
% ePair = [5, 7];
% elecLabels = {'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'};
% SESSION_ID = 'TremoLfpDBS-190927-100155';
nwbProcSuffix = ['_proc_lfp'];
nwbReadAcquPn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';
nwbWritePn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';

% tic
%%  Get bipolar lfp channel C0-C2 data for washin
clear raw epochIntervals

% read MASTER nwb file (with raw-acquisition data) 
nwbAcqu = nwbRead([nwbReadAcquPn SESSION_ID '_MASTER.nwb']);
raw = nwbAcqu.acquisition.get('rawVoltageTrace').deref();
fsRaw = raw.starting_time_rate;

decFactor = 24;
fsLfp = fsRaw / decFactor;


% % Collect sample lims for the multiple tdt "files" represented in NWB
% % all values to be found in the string within raw.description
% sampLims = [1, 414740480; ...
%             414740481, 502513664; ...
%             502513665, 590540800; ...
%             590540801, 676052992; ...
%             676052993, 764444672; ...
%             764444673, 811200512];
        
        
        
% Loop thru decimating raw ephys and filtering to LFP, storing in LFP
% matrix
lfpMat = [];
tstMat = [];
nFiles = length(sampLims);
nChans = 8;
% tic
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
% toc



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
% toc

% end



end



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