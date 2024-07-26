% script for converting acquisition accel data into NWB processing format. Done so that
% data between Mela, Kiwi, and Uva may be easily compared and the same
% scripts run on them all, since original formats are all frikkiin
% different. 

% 2 mg/kg day
% SESSION_ID = 'harPilotMela-150928-094900';

% 6 mg/kg day
% SESSION_ID = 'harPilotMela-150930-082000';

% 12 mg/kg day
% SESSION_ID = 'harPilotMela-151013-072100';



%% 

% Collected on microstrain G-link hardware
clear; close all

% % G-link sampling freq when streaming 3 channels (hardware design)
% fs = 617; % samples/second


% MANUAL INPUTS:
% projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';
% dataAcqPn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\';
% ACCdataPn = [projRootPath 'Data Acquisition\Kiwi Experiment\'];
% SESSION_ID_prefix = 'harPilotKiwi-';


% Fullpath for directory in which to save NWB file results
% nwbSavePn = [dataAcqPn 'NWBdata\KiwiPilot\'];
nwbSuffix = '_ACCproc';

nwbDataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';
nwbFullPath = [nwbDataPath 'MelaPilot\'];
projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';


SESSION_ID = 'harPilotMela-151013-072100';


% Load nwb master file for this session, get acc TimeSeries object
disp('Reading NWB file...')
nwb = nwbRead([nwbFullPath SESSION_ID '_ACCdata.nwb']);
disp('DONE!')

% 
% % % Speicfy which day's session-data is being concatenated and stored in NWB
% % % format
% % sessionNum = 3;
% 
% 
% % INPUT METADATA TABLES
% % Read in metadata table for Mela's recordings
% 
% tablePath = 'Data Acquisition\Kiwi Experiment\';
% metaTab = readtable([projRootPath tablePath 'RecordingsMetadata2.xlsx']);
% % harTab = readtable([projRootPath tablePath 'harmalineSessionData.csv']);
% fullMeta = metaTab;
% 
% 
% % Remove rows not pertaining to the experiment
% isBad = strcmp(fullMeta.trialType, 'test');
% fullMeta = fullMeta(~isBad,:);

%% Perform filtering and downsampling to get all Kiwi data on the same consistent sampling frequency

% Load all acc data 
accObj = nwb.acquisition.get('ACC');
acc = accObj.data.load;
tst = accObj.timestamps.load;
ctrl = accObj.control.load;



%% Create NWB object and Add in data to NWB object and save final result
% This file to be accessed by the "external link" property of the MASTER
% nwb file

% give file a globally unique id, more for machine readability than human
id = genUUID;

% Init NWB
nwbProc = NwbFile(...
    'session_start_time', nwb.session_start_time, ...
    'identifier', id, ...
    'session_description', 'a harmaline-recording pilot study', ...
    'timestamps_reference_time', nwb.timestamps_reference_time);

% Fill in the values to be included in the NWB file
time_series2 = types.core.TimeSeries( ...
    'comments', 'Data stored as double-precision columns of IMU readings over time.', ...
    'data', acc, ...
    'data_resolution', -1.0, ... % default for unknown
    'data_unit', 'missing', ...
    'description', 'AccX, AccY, AccZ', ...
    'timestamps', tst, ...
    'starting_time', 0, ...
    'starting_time_rate', 617, ...
    'control', ctrl, ...
    'control_description', 'numbers label the data samples correspond to N different recordings concatenated within this NWB file');

% nwbProc.acquisition.set('ACCproc', time_series2);

acc_module = types.core.ProcessingModule(...
    'description', 'holds accelerometer data, resampled to consistent sampling freq if needed.');
acc_module.nwbdatainterface.set('ACCproc', time_series2);
nwbProc.processing.set('ACCproc', acc_module);




% Write the nwb file in one go
tic

fullSavePath = [nwbFullPath SESSION_ID nwbSuffix '.nwb'];

% Next two lines essentially let you overwrite an existing version of this
% nwb file
nwbnrtl.io.checkdelete(fullSavePath); 
nwbExport(nwbProc, fullSavePath)
% nwbExport(nwb, [NWBPATH 'TEST' '.nwb'])

disp('DONE EXPORTING NWB!');
toc



%% Test loading and navigating the new nwb file

nwb2 = nwbRead(fullSavePath)



%% SUB-FUNCTIONS

function  [accRes, tstRes] = resampleData(acc, tst, fsDesired)
% Resample accelerometry
cfg.upSampFactor = 50;

% Define threshold for timestamp spacing inside which we would have to
% replace interpolated values with NaNs (note, necessary for really big
% gaps)
% intervThresh = cfg.resample.gapThresh; % seconds
dt_tst = diff(tst);
% isGap = dt_tst > intervThresh;

% desired resampling frequency is a factor of the LFP sampling frequency:
% 2,034.5 Hz, divide by 20
% fsDesired = cfg.resample.fsDesired; 
fsUpsamp = fsDesired * cfg.upSampFactor;

% disp('Resampling non-uniform IMU data...')
[sigUniform, tUniform] = util.upsampUniform(acc, tst, fsUpsamp);
        
for i = 1:3
    dataDec(:,i) = decimate(sigUniform(:,i), cfg.upSampFactor, 'fir');
    
end
timeDec = downsample(tUniform, cfg.upSampFactor);
disp('DONE!')

% % detrend the data of steady state and super slow movement
% fc = 0.1; % Hz highpass
% [b,a] = butter(2, fc / (fsDesired/2), 'high');
% dataFilt = filtfilt(b, a, dataDec);
% 
% accData = dataFilt;

accRes = dataDec;
tstRes = timeDec;

% % Assign NaN value to any data taking place within a long gap of time
% % according to original timestamps. Current data may be a simple linear
% % interpolation between a large time gap, so make it a NaN to keep analysis
% % good...
% nGaps = sum(isGap);
% idxGaps = find(isGap);
% isNanTime = false(length(accData), 1);
% for iGap = 1:nGaps
%     timeBand = [tst(idxGaps(iGap)), tst(idxGaps(iGap)+1)];
%     
%     % Discover which indices in newly-upsampled time and data pertain to gap times
%     temp = (timeDec >= timeBand(1)) & (timeDec < timeBand(2));
%     isNanTime = isNanTime | temp;
%     
% end
% % Set NaN values for accelerometry data interpolated between large gaps
% accData(isNanTime,1:3) = NaN;  
% 
% tstRes = timeDec;
% fs = fsDesired;
% 



end
function [isToChange] = identifysegtochange(accSeg, tstSeg, fsThresh)


nSegs = length(tstSeg);
isToChange = false(nSegs, 1);
for iSeg = 1:nSegs
    dt = median(diff(tstSeg{iSeg}));
    fs = 1/dt;
    
    % change those segments
    % Because of stupidity I recorded things at two different sampling
    % frequencies. Now, PSD peak amplitudes are significantly affected by 
    % sampling frequency, so to properly compare data recorded with
    % different fs we need to downsample the faster one to get as close to
    % the slower as possible. 

    % two frequencies used: 500 Hz, 8,193.4 Hz
    if fs > fsThresh
        isToChange(iSeg) = true;
        
    end
%     %     DWNFACTOR = 16;
%         fs = fs/cfg.DWNFACTOR; % bring down to 512.0852
%     %     t = downsample(t, cfg.DWNFACTOR);
%         clear y
%         for i = 1:3
%             y(:,i) = decimate(accDetrend(:,i), cfg.DWNFACTOR);
% 
%         end
%         accDetrend = y;
% 
%     else
%     %     dTime = t;
%     end
    
    



    
    
end

end

function [accSeg, tstSeg] = segmentacc(acc, tst, ctrl, maxSamp)
% separate acc data and timestamps into time segments based on recordings.
% Large recordings are further subdivided into smaller time segements according
% to max segment time. maxTime (seconds). ctrl has zero at the start of
% every recording. 

if ~exist('maxSamp'), maxSamp = []; end


% fs = 617;
% maxSamp = round(maxTime * fs); 
begIdx = find(ctrl == 0);
endIdx = [(begIdx(2:end)-1); length(ctrl)];
% nSamps = diff([begIdx endIdx], [], 2);

% enforce maximum samples for each segment
% isToobig = nSamps > maxSamp;
if isempty(maxSamp)
    begIdxCorr = begIdx;
    endIdxCorr = endIdx;
    
else
    nSamps = diff([begIdx endIdx], [], 2);
    begIdxCorr = [];
    endIdxCorr = [];
    for i = 1:length(begIdx)

        if diff([begIdx(i) endIdx(i)], [], 2) > maxSamp % Too big! break up the segment into smaller ones
            nRowsAdd = floor(nSamps(i)/maxSamp); 
            xb = [1:nRowsAdd]';

            % begIdxCorr
            begIdxSmall = begIdx(i) + xb*maxSamp + 1;
            begIdxCorr = [begIdxCorr; begIdx(i); begIdxSmall];

            % endIdxCorr
            endIdxSmall = begIdx(i) + xb*maxSamp;
            endIdxCorr = [endIdxCorr; endIdxSmall; endIdx(i)];

        else
            begIdxCorr = [begIdxCorr; begIdx(i)];
            endIdxCorr = [endIdxCorr; endIdx(i)];

        end



    end
    
end


% Gather data segments
nSegs = length(begIdxCorr);
for iSeg = 1:nSegs
    accSeg{iSeg,1} = acc(begIdxCorr(iSeg):endIdxCorr(iSeg),:);
    tstSeg{iSeg,1} = tst(begIdxCorr(iSeg):endIdxCorr(iSeg),:);
    
end


end

function [i_acc,i_tst,i_acq] = getAccInfoKiwi(data, time)
% if it's a stim recording, the first channel will be a recording of TACS
% or DBS signal, for a total of 4 channels; otherwise there should just be
% 3 channels of accel
if size(data, 2) == 4 
    i_acc = data(:,2:4);

else
    i_acc = data(:,1:3);

end
i_tst = time;
sampNum = length(i_tst);
i_acq = [0:(sampNum-1)]';
   
end

function descrString = buildElectricalSeriesDescr(TDTblock_multiple, TDTblock_sampsCumul)
descrString = 'Concatenated TDTblocks and their sample-indices:';

nBlocks = length(TDTblock_multiple);
for iBlock = 1:nBlocks
    descrString = [descrString, ' ', num2str(iBlock), ') ', ...
        TDTblock_multiple{iBlock}, '[', ...
        num2str(TDTblock_sampsCumul(iBlock,1)), ',',...
        num2str(TDTblock_sampsCumul(iBlock,2)), '];']
    
end

end

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


