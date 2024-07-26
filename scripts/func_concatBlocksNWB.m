function nwbSess = func_concatBlocksNWB_export(nwbBlk, nwbexportPath)
% Function is specific to this project's data pipeline. 
% Takes in an array of NWB objects that each represent a conversion of TDT
% block data into NWB format, stored in "acquisistion". Objects also have the results of DBS artifact
% products of DBS artifact rejection and spike filtering steps, stored in
% "processing". 
% Concatenates all data with time referenced to the start time of the first
% TDT block in the sequence, creating a single file for all data in the
% recording session (Session contains 12 DBS epochs, different DBS
% conditions).
%
% ASSUMPTIONS: 
% Code assumes that identical fields are recorded (and exist) across
% all NWB objects to be concatenated.
% Code also assumes that the etirety of data in stores across all NWB objects
% to be concatenated can be held in RAM in one matlab session (Could be a
% problem for many big TDT blocks...).
%
% NOTE: Currently only "stream" type stores are translated by this function. All
% such stores are organized under the "acquisition" segment of the NWB
% object, with original store labels as detected by the Matlab SDK provided
% by TDT (TDTbin2mat.m). All data are stored as TimeSeries objects.
%
% INPUT
% TDTblocks - Cell array of strings, each cell containing a string of the
% full path to each TDTblock. Arranged in the order to be concatenated. 
%
% OUTPUT
% nwb - matlab NWB object with all data from TDT block
%
% NOTE on the output: any concatenated timeseries will have an empty field
% where "starting_time" used to be. This is to prevent confusion about
% whether the concatenated data is continuous with no breaks (which it is
% likely not), or whether the data is broken up by gaps between each block
% (which it most likely is). The "timestamps" field will show accurate
% times despite jumps in the data.


tic

% % Create an NWB object for each TDT block to be concatenated
% nBlocks = length(TDTblocks);
% for iBlk = 1:nBlocks
%     nwbBlk(iBlk) = contrib.tdt.TDTbin2NWB(TDTblocks{iBlk});
% 
% end



%% Concatenate all nwbBlk objects into one nwb object:
% nwbBlk must be an array of nwb objects, in order of occurrence
% Takes in the array of individual nwb objects and concatenates all data in
% the acquisition field.

% Initialize the session-wide NWB object
% info based on first nwb object in sequence  

nwbSess = NwbFile(...
    'session_start_time', nwbBlk(1).session_start_time, ...
    'general_session_id', nwbBlk(1).general_session_id, ...
    'session_description', nwbBlk(1).session_description, ...
    'identifier', nwbBlk(1).identifier, ...
    'timestamps_reference_time', nwbBlk(1).session_start_time);

% Export "first draft" of concatenated NWB file
[nwbSfilepath, nwbSname, ~] = fileparts(nwbexportPath);
nwbexportPath = [nwbSfilepath '\' nwbSname '.nwb'];
nwbExport(nwbSess, nwbexportPath);



%% Run concatenation process for each TimeSeries in acquisition

% Identify all TimeSeries objects to be concatenated
stores = nwbBlk(1).acquisition.keys;
nStores = length(stores);

nBlocks = length(nwbBlk);

% FOR EACH TIMESERIES
for iStore = 1:nStores
    disp(['Concatenating TimeSeries ', num2str(iStore), ...
        ' of ', num2str(nStores), ' (', stores{iStore}, ')...']);

    % Get single concatenated TimeSeries
%     nwbTSPath = join(["acquisition.get('" stores{iStore} "')"], "");
    for iBlk = 1:nBlocks
        tSeriesBlk(iBlk) = nwbBlk(iBlk).acquisition.get(stores{iStore});
        session_start_time(iBlk) = nwbBlk(iBlk).session_start_time;
        
    end
    tSeriesConc = util.concat_TimeSeries(tSeriesBlk, session_start_time);

    % Set the concat'd TimeSeries in the concat'd nwb file, and overwrite 
    nwbSess.acquisition.set(stores{iStore}, tSeriesConc);
    nwbSess = util.nwbOverwrite(nwbSess, nwbexportPath);
    
   disp(['TimeSeries ' stores{iStore} ' updated']);

end



%% Run concatenation process for each TimeSeries in acquisition
% 
% % Identify all TimeSeries objects to be concatenated
% stores = nwbBlk(1).acquisition.keys;
% nStores = length(stores);
% 
% nBlocks = length(nwb
% 
% % FOR EACH TIMESERIES
% for iStore = 1:nStores
%     disp(['Concatenating TimeSeries ', num2str(iStore), ...
%         ' of ', num2str(nStores), ' (', stores{iStore}, ')...']);
% 
%     % Get single concatenated TimeSeries
% %     nwbTSPath = join(["acquisition.get('" stores{iStore} "')"], "");
%     for 
%     
%     tSeriesConc = util.concat_TimeSeries(tSeriesBlk, session_start_time);
% 
%     % Set the concat'd TimeSeries in the concat'd nwb file, and overwrite 
%     nwbSess.acquisition.set(stores{iStore}, tSeriesConc);
%     nwbSess = util.nwbOverwrite(nwbSess, nwbexportPath);
%     
%    disp(['TimeSeries ' stores{iStore} ' updated']);
% 
% end



%% Run concatenation process for each Timeseries in acquisition (i.e. each TDT stream store)
% NOTE: Code is set up to clear data once it's been concatenated, in order
% to save RAM


% Preserve non-concatenating info for final concat'd TimeSeries, based on
% first block, for each store:
for iStore = 1:nStores
    fsSess(iStore) = nwbBlk(1).acquisition.get(stores(iStore)).starting_time_rate;
    srtTimeSess(iStore) = nwbBlk(1).acquisition.get(stores(iStore)).starting_time;

end

% FOR EACH TimeSeries:
nBlocks = length(nwbBlk);
timeStamps = cell(1,nStores);
data =  cell(1,nStores);
refDT = datetime(nwbSess.session_start_time);

% For each TDTblock
tags = cell(nBlocks,1);

tsTracking = cell(nBlocks, 1);
idxStart = zeros(1, nStores);

for iBlk = 1:nBlocks
    disp(['Concatenating block ' num2str(iBlk) ' of ' num2str(nBlocks) '...']);
    
    % Obtain TDTblock identifier, start-time, and stop-time for each TDT block 
    tags{iBlk,1} = nwbBlk(iBlk).general_session_id;
    
    refSec = seconds(datetime(nwbBlk(iBlk).timestamps_reference_time) - refDT);
    start_time(iBlk,1) = refSec;
    durstr = nwbBlk(iBlk).scratch.get('session_duration').data.load;
    durSec = seconds(datetime(durstr) - datetime('00:00:00'));
    stop_time(iBlk,1) = start_time(iBlk,1) + durSec; 
    
    
    % Iteratively concatenate timestamps and data
    for iStore = 1:nStores

        storeStr = stores{iStore};

        tSeries = nwbBlk(iBlk).acquisition.get(storeStr);

        % store the data
        startTime = tSeries.starting_time;
        fs = tSeries.starting_time_rate;
        tst = ((0:(size(tSeries.data,1)-1)) * (1/fs)) + startTime;


        % Re-reference block's timestamps relative to session timestamp reference time
        tstSync = tst + refSec;


        % concatenate timestamps
        timeStamps{iStore} = [timeStamps{iStore}; tstSync'];

        % concatenate data
        data{iStore} = [data{iStore}; tSeries.data];


        % Remove the timeseries object from current nwb block object, to save
        % memory
%         nwbBlk(iBlk).acquisition.remove(storeStr);
        clear tSeries
        
        % track timeseries data
        nSamps = length(tstSync);
        tss(iStore).idx_start = idxStart(iStore);
        tss(iStore).count = nSamps;
        tss(iStore).timeseries = types.untyped.ObjectView( ...
           ['/acquisition/' storeStr]);
       
        
        idxStart(iStore) = idxStart(iStore) + nSamps; 
       
    end
    
    tsTracking{iBlk,1} = tss;


end % END LOOP

% store the data as final TimeSeries object in the session nwb object
% Store data as an ElectricalSeries object
for iStore = 1:nStores
    
    % Save concatenated data into NWB session object
    tSeriesName = stores(iStore);
    timeSeries_sess = types.core.TimeSeries( ...
        'data', data{iStore}, ...
        'data_resolution', -1.0, ... % default for unknown
        'data_unit', '<missing>', ...
        'starting_time_rate', fsSess(iStore), ...
        'timestamps', timeStamps{iStore});

    nwbSess.acquisition.set(tSeriesName, timeSeries_sess);
    
    % Now that data was saved, clear duplicate data from RAM
    data{iStore} = [];
    timeStamps{iStore} = [];

end

% clear nwbBlk

clear fsSess srtTimeSess start_time stop_time data timeStamps







%% Run concatenation process for each Timeseries in acquisition (i.e. each TDT stream store)
% NOTE: Code is set up to clear data once it's been concatenated, in order
% to save RAM

nStores = length(stores);

% Preserve non-concatenating info for final concat'd TimeSeries, based on
% first block, for each store:
for iStore = 1:nStores
    fsSess(iStore) = nwbBlk(1).acquisition.get(stores(iStore)).starting_time_rate;
    srtTimeSess(iStore) = nwbBlk(1).acquisition.get(stores(iStore)).starting_time;

end

% FOR EACH TimeSeries:
nBlocks = length(nwbBlk);
timeStamps = cell(1,nStores);
data =  cell(1,nStores);
refDT = datetime(nwbSess.session_start_time);

% For each TDTblock
tags = cell(nBlocks,1);

tsTracking = cell(nBlocks, 1);
idxStart = zeros(1, nStores);

for iBlk = 1:nBlocks
    disp(['Concatenating block ' num2str(iBlk) ' of ' num2str(nBlocks) '...']);
    
    % Obtain TDTblock identifier, start-time, and stop-time for each TDT block 
    tags{iBlk,1} = nwbBlk(iBlk).general_session_id;
    
    refSec = seconds(datetime(nwbBlk(iBlk).timestamps_reference_time) - refDT);
    start_time(iBlk,1) = refSec;
    durstr = nwbBlk(iBlk).scratch.get('session_duration').data.load;
    durSec = seconds(datetime(durstr) - datetime('00:00:00'));
    stop_time(iBlk,1) = start_time(iBlk,1) + durSec; 
    
    
    % Iteratively concatenate timestamps and data
    for iStore = 1:nStores

        storeStr = stores{iStore};

        tSeries = nwbBlk(iBlk).acquisition.get(storeStr);

        % store the data
        startTime = tSeries.starting_time;
        fs = tSeries.starting_time_rate;
        tst = ((0:(size(tSeries.data,1)-1)) * (1/fs)) + startTime;


        % Re-reference block's timestamps relative to session timestamp reference time
        tstSync = tst + refSec;


        % concatenate timestamps
        timeStamps{iStore} = [timeStamps{iStore}; tstSync'];

        % concatenate data
        data{iStore} = [data{iStore}; tSeries.data];


        % Remove the timeseries object from current nwb block object, to save
        % memory
%         nwbBlk(iBlk).acquisition.remove(storeStr);
        clear tSeries
        
        % track timeseries data
        nSamps = length(tstSync);
        tss(iStore).idx_start = idxStart(iStore);
        tss(iStore).count = nSamps;
        tss(iStore).timeseries = types.untyped.ObjectView( ...
           ['/acquisition/' storeStr]);
       
        
        idxStart(iStore) = idxStart(iStore) + nSamps; 
       
    end
    
    tsTracking{iBlk,1} = tss;


end % END LOOP

% store the data as final TimeSeries object in the session nwb object
% Store data as an ElectricalSeries object
for iStore = 1:nStores
    
    % Save concatenated data into NWB session object
    tSeriesName = stores(iStore);
    timeSeries_sess = types.core.TimeSeries( ...
        'data', data{iStore}, ...
        'data_resolution', -1.0, ... % default for unknown
        'data_unit', '<missing>', ...
        'starting_time_rate', fsSess(iStore), ...
        'timestamps', timeStamps{iStore});

    nwbSess.acquisition.set(tSeriesName, timeSeries_sess);
    
    % Now that data was saved, clear duplicate data from RAM
    data{iStore} = [];
    timeStamps{iStore} = [];

end

% clear nwbBlk

clear fsSess srtTimeSess start_time stop_time data timeStamps


% %% Run concatenation process for each Timeseries in PROCESSING
% % NOTE: Code is set up to clear data once it's been concatenated, in order
% % to save RAM
% 
% % Identify all TimeSeries objects to be concatenated
% procMod = 'DBSartifact_detect-clean';
% % stores = nwbBlk(1).processing.get(procMod).nwbdatainterface.keys;
% stores = 'spkClean';
% 
% % nStores = length(stores);
% nStores = 1;
% 
% % Preserve non-concatenating info for final concat'd TimeSeries, based on
% % first block, for each store:
% for iStore = 1:nStores
%     fsSess(iStore) = nwbBlk(1).processing.get(procMod).nwbdatainterface.get(stores).starting_time_rate;
%     srtTimeSess(iStore) = nwbBlk(1).processing.get(procMod).nwbdatainterface.get(stores).starting_time;
% 
% end
% 
% % FOR EACH TimeSeries:
% nBlocks = length(nwbBlk);
% timeStamps = cell(1,nStores);
% data =  cell(1,nStores);
% refDT = datetime(nwbSess.session_start_time);
% 
% % For each TDTblock
% tags = cell(nBlocks,1);
% 
% tsTracking = cell(nBlocks, 1);
% idxStart = zeros(1, nStores);
% 
% for iBlk = 1:nBlocks
%     disp(['Concatenating block ' num2str(iBlk) ' of ' num2str(nBlocks) '...']);
%     
%     % Obtain TDTblock identifier, start-time, and stop-time for each TDT block 
%     tags{iBlk,1} = nwbBlk(iBlk).general_session_id;
%     
%     refSec = seconds(datetime(nwbBlk(iBlk).timestamps_reference_time) - refDT);
%     start_time(iBlk,1) = refSec;
%     durstr = nwbBlk(iBlk).scratch.get('session_duration').data.load;
%     durSec = seconds(datetime(durstr) - datetime('00:00:00'));
%     stop_time(iBlk,1) = start_time(iBlk,1) + durSec; 
%     
%     
%     % Iteratively concatenate timestamps and data
%     for iStore = 1:nStores
% 
%         storeStr = stores;
% 
%         tSeries = nwbBlk(iBlk).processing.get(procMod).nwbdatainterface.get(stores);
% 
%         % store the data
%         startTime = tSeries.starting_time;
%         fs = tSeries.starting_time_rate;
%         tst = ((0:(size(tSeries.data,1)-1)) * (1/fs)) + startTime;
% 
% 
%         % Re-reference block's timestamps relative to session timestamp reference time
%         tstSync = tst + refSec;
% 
% 
%         % concatenate timestamps
%         timeStamps{iStore} = [timeStamps{iStore}; tstSync'];
% 
%         % concatenate data
%         data{iStore} = [data{iStore}; tSeries.data];
% 
% 
%         % Remove the timeseries object from current nwb block object, to save
%         % memory
% %         nwbBlk(iBlk).acquisition.remove(storeStr);
%         clear tSeries
%         
% %         % track timeseries data
% %         nSamps = length(tstSync);
% %         tss(iStore).idx_start = idxStart(iStore);
% %         tss(iStore).count = nSamps;
% %         tss(iStore).timeseries = types.untyped.ObjectView( ...
% %            ['/acquisition/' storeStr]);
%        
%         
% %         idxStart(iStore) = idxStart(iStore) + nSamps; 
%        
%     end
%     
% %     tsTracking{iBlk,1} = tss;
% 
% 
% end % END LOOP
% 
% % store the data as final TimeSeries object in the session nwb object
% % Store data as an ElectricalSeries object
% 
% ecephys_module = types.core.ProcessingModule(...
%     'description', 'Holds ephys data that has been DBS artifact-cleaned');
% % ecephys_module.nwbdatainterface.set('tDBS_uncorrected', tDBS_uncorrected);
% % ecephys_module.nwbdatainterface.set('tDBS', tDBS);
% % ecephys_module.nwbdatainterface.set('tDBS_virtPre', tDBS_virtPre);
% % ecephys_module.nwbdatainterface.set('tDBS_virtPos', tDBS_virtPos);
% % ecephys_module.nwbdatainterface.set('spkClean', TSspkSclean);
% 
% 
% 
% for iStore = 1:nStores
%     
%     % Save concatenated data into NWB session object
%     tSeriesName = stores;
%     timeSeries_sess = nwbBlk(iBlk).processing.get(procMod).nwbdatainterface.get(tSeriesName);
%     timeSeries_sess.data = data{iStore};
%     timeSeries_sess.timestamps = timeStamps{iStore};
% %     timeSeries_sess = types.core.TimeSeries( ...
% %         'data', data{iStore}, ...
% %         'data_resolution', -1.0, ... % default for unknown
% %         'data_unit', '<missing>', ...
% %         'starting_time_rate', fsSess(iStore), ...
% %         'timestamps', timeStamps{iStore});
% 
%     ecephys_module.nwbdatainterface.set(tSeriesName, timeSeries_sess);
% 
% 
%     nwbSess.acquisition.set(tSeriesName, timeSeries_sess);
%     nwbSess.processing.set(procMod).nwbdatainterface.set(stores)
%     % Now that data was saved, clear duplicate data from RAM
%     data{iStore} = [];
%     timeStamps{iStore} = [];
% 
% end
% 
% % clear nwbBlk
% nwb.processing.set('DBSartifact_detect-clean', ecephys_module);



%% Create TimeIntervals object to partition data by TDT block
% Allows for fast loading access if all data from one block is desired.

epochs = types.core.TimeIntervals( ...
    'id', types.hdmf_common.ElementIdentifiers( ...
        'data', 0:(nBlocks-1)), ...
    'colnames', {'tags', 'start_time','stop_time'}, ...
    'description', 'time intervals for original epoch partitions during data acquisition.', ...
    'start_time', types.hdmf_common.VectorData( ...
        'data', start_time, ...
        'description','start time of epoch'), ...
    'stop_time', types.hdmf_common.VectorData( ...
        'data', stop_time, ...
        'description','end time of epoch'), ...
    'tags', types.hdmf_common.VectorData( ...
        'data', tags, ...
        'description', 'Labels of separate experimental epochs'));
    
nwbSess.intervals_epochs = epochs;

disp('CONCATENATION COMPLETE!');
toc
% BELOW IS THE ORIGINAL ATTEMPT TO INCORPORATE TIMESERIES OBJECT REFERENCES
% Unsure yet how to set this up properly, initial attempts have failed...
% epochs = types.core.TimeIntervals( ...
%     'id', types.hdmf_common.ElementIdentifiers( ...
%         'data', 0:(nBlocks-1)), ...
%     'colnames', {'tags', 'start_time','stop_time', 'timeseries'}, ...
%     'description', 'time intervals for original epoch partitions during data acquisition.', ...
%     'start_time', types.hdmf_common.VectorData( ...
%         'data', start_time, ...
%         'description','start time of epoch'), ...
%     'stop_time', types.hdmf_common.VectorData( ...
%         'data', stop_time, ...
%         'description','end time of epoch'), ...
%     'tags', types.hdmf_common.VectorData( ...
%         'data', tags, ...
%         'description', 'Labels of separate experimental epochs, may not be in chron order.'), ...
%     'timeseries', types.hdmf_common.VectorData( ...
%         'data', tsTracking, ...
%         'description', 'timeseries stuff'));
%     
% nwbSess.intervals_epochs = epochs;



end

% SUB-FUNCTIONS

function nwbSess = nwbconcat(nwbBlk)
% nwbBlk must be an array of nwb objects, in order of occurrence
% Takes in the array of individual nwb objects and concatenates all data in
% the acquisition field.

% Identify all TimeSeries objects to be concatenated
stores = nwbBlk(1).acquisition.keys;



% Initialize the session-wide NWB object
% info based on first nwb object in sequence  
nwbSess = NwbFile(...
    'session_start_time', nwbBlk(1).session_start_time, ...
    'identifier', nwbBlk(1).identifier, ...
    'session_description', nwbBlk(1).session_description, ...
    'timestamps_reference_time', nwbBlk(1).session_start_time);



%% Run concatenation process for each Timeseries in acquisition (i.e. each TDT stream store)
% NOTE: Code is set up to clear data once it's been concatenated, in order
% to save RAM



% FOR EACH TimeSeries:
nStores = length(stores);
for iStore = 1:nStores
    storeStr = stores{iStore};

    % Preserve non-concatenating info for final concat'd TimeSeries, based on
    % first block:
    fsSess = nwbBlk(1).acquisition.get(storeStr).starting_time_rate;
    srtTimeSess = nwbBlk(1).acquisition.get(storeStr).starting_time;


    % Iteratively concatenate timestamps and data
    timeStamps = [];
    data = [];
    refDT = datetime(nwbSess.session_start_time);
    nBlocks = length(nwbBlk);
    for iBlk = 1:nBlocks

        tSeries = nwbBlk(iBlk).acquisition.get(storeStr);

        % store the data
        startTime = tSeries.starting_time;
        fs = tSeries.starting_time_rate;
        tst = ((0:(size(tSeries.data,1)-1)) * (1/fs)) + startTime;


        % Re-reference block's timestamps relative to session timestamp reference time
        refSec = seconds(datetime(nwbBlk(iBlk).timestamps_reference_time) - refDT);
        tstSync = tst + refSec;


        % concatenate timestamps
        timeStamps = [timeStamps; tstSync'];

        % concatenate data
        data = [data; tSeries.data];


        % Remove the timeseries object from current nwb block object, to save
        % memory
        nwbBlk(iBlk).acquisition.remove(storeStr);
        clear tSeries

    end

    % store the data as final TimeSeries object in the session nwb object
    % Store data as an ElectricalSeries object
    tSeriesName = storeStr;
    timeSeries_sess = types.core.TimeSeries( ...
        'data', data, ...
        'data_resolution', -1.0, ... % default for unknown
        'starting_time', srtTimeSess, ... % seconds
        'starting_time_rate', fsSess, ...
        'timestamps', timeStamps);

    nwbSess.acquisition.set(tSeriesName, timeSeries_sess);

end

end

