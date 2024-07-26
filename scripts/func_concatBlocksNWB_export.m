function nwbSess = func_concatBlocksNWB_export(nwbBlk, nwbexportPath)
% Function is specific to this project's data pipeline. 
% Takes in an array of NWB objects that each represent a conversion of TDT
% block data into NWB format. Objects also have the results of DBS artifact
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
% Note: the util.nwbOverwrite code allows for dealing with only one
% TimeSeries concatenation at a time, so as not to overwhelm RAM. 
%
%
% INPUT
% nwbBlk - Array of NwbFile objects, arranged in the order to be concatenated. 
% nwbexportPath - full path to where the concat'd NWB object is to be
% saved. Including the name for the file and .nwb extension.
%
% OUTPUT
% nwbSess - matlab NWB object with all concatenated data. 
%
% NOTE on the output: As concatenated data is likely to have empty time
% skips between TDTblock recordings, TimeSeries will have timestamps for
% reconstructing the time sample of each value relative to start of
% session. This is the correct way to track time in concatenated data like
% this. Don't use the sampling rate and starting time to get time data! 

tic

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


% special section for tracking samples of spkClean
nBlocks = length(nwbBlk);
idx_start = zeros(nBlocks, 1);
 idx_stop = zeros(nBlocks, 1);
sampsCum = 0;
procMod = 'DBSartifact_detect-clean';
for iBlk = 1:nBlocks
    tSeriesBlkSpk(iBlk) = nwbBlk(1).processing.get(procMod).nwbdatainterface.get('spkClean');
    if isa(tSeriesBlkSpk(iBlk).data, 'types.untyped.DataStub')
        N = tSeriesBlkSpk(iBlk).data.dims(1);
        
    elseif isnumeric(tSeriesBlkSpk(iBlk).data)
        N = size(tSeriesBlkSpk(iBlk).data, 1);
        
    else
        error('check this out')
        
    end
    idx_start(iBlk,1) = sampsCum + 1;
     idx_stop(iBlk,1) = sampsCum + N;
    
    sampsCum = sampsCum + N;
        
end


%% Run concatenation process for each TimeSeries in acquisition
disp('For NWB ACQUISITION section: ');

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
    
   disp(['TimeSeries updated: ' stores{iStore}]);

end

disp('ACQUISITION section DONE');



%% Run concatenation process for each TimeSeries in acquisition
disp('For NWB PROCESSING section: ');

% Initialize processing module for concat version
ecephys_module = types.core.ProcessingModule(...
    'description', 'Holds ephys data that has been DBS artifact-cleaned');
procMod = 'DBSartifact_detect-clean';

% Identify all TimeSeries objects to be concatenated
stores = nwbBlk(1).processing.get(procMod).nwbdatainterface.keys;
nStores = length(stores);
nBlocks = length(nwbBlk);

% FOR EACH TIMESERIES
for iStore = 1:nStores
    disp(['Concatenating TimeSeries ', num2str(iStore), ...
        ' of ', num2str(nStores), ' (', stores{iStore}, ')...']);

    % Get single concatenated TimeSeries
%     nwbTSPath = join(["acquisition.get('" stores{iStore} "')"], "");
    for iBlk = 1:nBlocks
        tSeriesBlk(iBlk) = nwbBlk(iBlk).processing.get(procMod).nwbdatainterface.get(stores{iStore});
        session_start_time(iBlk) = nwbBlk(iBlk).session_start_time;
        
    end
    tSeriesConc = util.concat_TimeSeries(tSeriesBlk, session_start_time);

    % Set the concat'd TimeSeries in the concat'd nwb file, and overwrite 
    ecephys_module.nwbdatainterface.set(stores{iStore}, tSeriesConc);
    
   disp(['TimeSeries updated: ' stores{iStore}]);

end


nwbSess.processing.set(procMod, ecephys_module);
nwbSess = util.nwbOverwrite(nwbSess, nwbexportPath);


disp('PROCESSING section DONE');



%% Create TimeIntervals object to partition data by TDT block
% Allows for fast loading access if all data from one block is desired.
disp('Including times for each original block in intervals/epochs as a TimeIntervals object')
refDT = datetime(nwbSess.session_start_time);

nBlocks = length(nwbBlk);
for iBlk = 1:nBlocks
    refSec = seconds(datetime(nwbBlk(iBlk).timestamps_reference_time) - refDT);
    start_time(iBlk,1) = refSec;
    durstr = nwbBlk(iBlk).scratch.get('session_duration').data.load;
    durSec = seconds(datetime(durstr) - datetime('00:00:00'));
    stop_time(iBlk,1) = start_time(iBlk,1) + durSec; 
    tags{iBlk,1} = nwbBlk(iBlk).general_session_id;
    
end


    
% Create the TimeIntervals object for tracking original blocks intervals
epochs = types.core.TimeIntervals( ...
    'id', types.hdmf_common.ElementIdentifiers( ...
        'data', 0:(nBlocks-1)), ...
    'colnames', {'tags', 'start_time','stop_time', 'spkClean_idx_start', 'spkClean_idx_stop'}, ...
    'description', 'time intervals for original epoch partitions during data acquisition.', ...
    'start_time', types.hdmf_common.VectorData( ...
        'data', start_time, ...
        'description','start time of epoch'), ...
    'stop_time', types.hdmf_common.VectorData( ...
        'data', stop_time, ...
        'description','end time of epoch'), ...
    'spkClean_idx_start', types.hdmf_common.VectorData( ...
        'data', idx_start, ...
        'description', 'tracking the samples in spkClean for each epoch'), ...
    'spkClean_idx_stop', types.hdmf_common.VectorData( ...
        'data', idx_stop, ...
        'description', 'tracking the samples in spkClean for each epoch'), ...
    'tags', types.hdmf_common.VectorData( ...
        'data', tags, ...
        'description', 'Labels of separate experimental epochs'));
    
    
nwbSess.intervals_epochs = epochs;
nwbSess = util.nwbOverwrite(nwbSess, nwbexportPath);

disp('CONCATENATION & EXPORT COMPLETE!');
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