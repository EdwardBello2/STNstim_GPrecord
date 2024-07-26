% Convert all relevant data from TDT blocks from entire recording session
% into one nwb file
%
%

% Authro: Ed Bello
% Created: 3/26/2020
%
%
% ChangeLog:
% 
%
% To-Do
% - use UUID for globally unique file identifier
%
%
%% Call toolboxes

% remove TDTMatlabSDK matnwb version, and instead include my local version
% This is because matnwb current version has a bug when I try to read the
% harmaline tdt blocks, has to do with me having a TDT note file in the
% block

% This code assumes that NMRCNWBStandardization (and subdirectories) has 
% already been added to matlab path
% 
% addpath(genpath('C:\Users\bello043\Documents\GitHub\NMRCNWBStandardization'));
% 
% warning('off','all');
% rmpath(genpath('C:\Users\bello043\Documents\GitHub\NMRCNWBStandardization\toolbox\TDTMatlabSDK'));
% warning('on','all');
% addpath(genpath('C:\Users\bello043\Documents\GitHub\TDTMatlabSDK'));
% 
% 
% % Add project folder and all subdirectories within:
% addpath(genpath('C:\Users\bello043\Documents\GitHub\Thalamic-DBS-HarmalineTremor'));
% 
% % Remove .git portion of the file (some .git folders are huge/not needed
% % here), and suppress the MATLAB warning if the project folder is not
% % git-source-controlled:
% % warning('off', 'MATLAB:rmpath:DirNotFound');
% % rmpath([projRootPath, '\.git']);
% 
% % Add path to TDTmatlabSDK, which contains TDTbin2mat (loading TDT data)
% tdtmatlabsdkPath = 'C:\Users\bello043\Documents\GitHub\TDTMatlabSDK';
% addpath(genpath(tdtmatlabsdkPath));
% rmpath([tdtmatlabsdkPath, '\.git']);


% Add path to Git version-controlled code for Harmaline project
addpath(genpath('C:\Users\bello043\Documents\GitHub\Thalamic-DBS-HarmalineTremor'));



%% CONSTANTS
% Note: nwb file will be named identically to first TDT block of the
% session

% TDT related Data
DATALOCPATH = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\'; 
TDTTANKPATH = 'Uva-191115\'; 
FIRST_BLOCK = 'TremoLfpDBS-191115-100127';
TDTBLOCKS = {'TremoLfpDBS-191115-100127', 'TremoLfpDBS-191115-122420', ...
    'TremoLfpDBS-191115-125440', 'TremoLfpDBS-191115-132410', ...
    'TremoLfpDBS-191115-140755', 'TremoLfpDBS-191115-143421'};
NHP_RARID = '13LP1';

STORE = 'RAW8'; 
TDT_CH = [1:8];


% NWB save
SAVEPATH = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWB\';



%% Initialize nwb object based on TDT metadata and user input

fulltdtpath = [DATALOCPATH, TDTTANKPATH, FIRST_BLOCK];
tdt = TDTbin2mat(fulltdtpath, 'TYPE', {'epocs', 'snips', 'scalars'});
dt = datetime(tdt.info.date);
tt = datetime(tdt.info.utcStartTime); 
tott = dt + timeofday(tt);
sess_utcStart = datestr(tott, 'yyyy-mm-ddTHH:MM:SSZ');
% id = [sess_utcStart(1:4), sess_utcStart(6:7), sess_utcStart(9:10), '_', ...
%       NHP_RARID];

% give file a globally unique id, more for machine readability than human
id = genUUID;

  
nwb = nwbfile(...
    'session_start_time', sess_utcStart, ...
    'identifier', id, ...
    'session_description', 'a baseline day recording', ...
    'timestamps_reference_time', sess_utcStart);



% Fill "general" section of nwb file, constants set inside function 
nwb = nwbfill_general(nwb);

% Fill "electrode" information of nwb file, constants set inside function
nwb = nwbfill_electrode(nwb);


% Save initial nwb file so far
% ("lazy loading" property of nwb file means all data is not loaded at
% once, only that which is specifically called, saving RAM

% nwbExport(nwb, [SAVEPATH FIRST_BLOCK '.nwb'])
% clear nwb
% 
% % Now read it back in to test
% nwb = nwbRead([SAVEPATH FIRST_BLOCK '.nwb']);
% 
% 

%% Add raw TDT data in to acquisition

% % Recording of interest:
% FULL = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata';
% TDTTANKPATH = 'Uva-191115';
% TDTBLOCKS = {'TremoLfpDBS-191115-100127', 'TremoLfpDBS-191115-122420', ...
%     'TremoLfpDBS-191115-125440', 'TremoLfpDBS-191115-132410', ...
%     'TremoLfpDBS-191115-140755', 'TremoLfpDBS-191115-143421'};
nBlks = length(TDTBLOCKS);
fullTankPath = [DATALOCPATH TDTTANKPATH(1:end-1)]; % remove the "\" from the end of string
  
for ch = TDT_CH
%     % Test case:
%     tdtConcat.streams.(STORE).data = 1:1000;
%     fsData = 1000;
%     t_array = 1:1000;
    
    
    % Concatenate data from all TDTblocks in the session, and get time
    % samples to match data samples with proper discontinuities
    [tdtConcat, idxBlk, idxBlkSync, infoBlk] = concatTDTstream_singlechan(fullTankPath, TDTBLOCKS, STORE, ch);
    fsData = tdtConcat.streams.(STORE).fs;
    
    t_array = zeros(idxBlk(end), 1);
    nBlks = size(idxBlk, 1);
    for iBlk = 1:nBlks
            t_array(idxBlk(iBlk,1):idxBlk(iBlk,2)) = (1/fsData) * ((idxBlkSync(iBlk,1):idxBlkSync(iBlk,2)) - 1);

    end


    % specify the electrode(s) that this signal came from according to nwb
    electrodes_object_view = types.untyped.ObjectView( ...
        '/general/extracellular_ephys/electrodes');

    electrode_table_region = types.core.DynamicTableRegion( ...
        'table', electrodes_object_view, ...
        'description', ['C' num2str(ch-1)], ...
        'data', (ch - 1)); % nwb uses zero as first index...


    % store it as electrical series object in nwb
    electrical_series = types.core.ElectricalSeries( ...
        'comments', 'data stored as single-precision matlab array, indices for mapping current data to previous TDT blocks is in "description"', ...
        'data', tdtConcat.streams.(STORE).data, ...
        'data_resolution', -1.0, ... % default for unknown
        'data_unit', 'Volts', ...
        'description', ['ch' num2str(ch) '-C' num2str(ch-1)], ...
        'electrodes', electrode_table_region, ... % points to value in electrode table
        'starting_time', 0.0, ... % seconds
        'starting_time_rate', fsData);
%         'starting_time_unit', 'Seconds');
%         'timestamps', t_array, ...
%         'timestamps_interval', 1, ...
%         'timestamps_unit', 'Seconds');

%     nwb.acquisition.set(['ElectricalSeries_ch', num2str(ch)], electrical_series);
    ElSerStr = ['ElectricalSeries_ch' num2str(ch)];
    nwb.acquisition.set(ElSerStr, electrical_series);

%     
%     % to save RAM, export current nwb, clear it, then load nwb again
%     nwbExport(nwb, [SAVEPATH FIRST_BLOCK '.nwb'])
%     clear nwb
%     nwb = nwbRead([SAVEPATH FIRST_BLOCK '.nwb']);
    
end

nwbExport(nwb, [SAVEPATH FIRST_BLOCK '.nwb'])

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
    'date_of_birth', '2001', 'description', 'nickname: Uva, source: breeder', ...
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
