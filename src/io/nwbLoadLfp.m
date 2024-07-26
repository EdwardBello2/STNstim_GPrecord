function [lfp, timestamps, fs] = nwbLoadLfp(nwb, varargin)
% default is just to load all data at once

% Load only data within a certain epoch

defaultEpoch = 'all';
defaultChan = [];


p = inputParser;

addRequired(p, 'nwb');
addParameter(p, 'epoch', defaultEpoch, @ischar);
addParameter(p, 'channel', defaultChan, @isnumeric);
parse(p, nwb, varargin{:});


elecSeries = nwb.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('lfp');

nChans = elecSeries.data.dims(2);
nSamps = elecSeries.data.dims(1);



ch = p.Results.channel;
epoch = p.Results.epoch;
nwb = p.Results.nwb;



%% Determine sample intervals to load based on epoch specified
% if not specified, load all

% Specify range of samples to load from LFPs based on epoch
switch epoch
    case 'all'
        samps = [1, nSamps];
        
    case {'naiveWashinBase', 'naiveTaskBase', 'naiveDbsBase', ...
            'naivePreExp', 'naiveDbsExp', 'harWashinExp', ...
            'harTaskExp', 'harDbsExp','harWashoutExp'}
        % Get the time-interval of this epoch, in secods
        tbl = nwb2table_epochs(nwb); % times in this table are in seconds, by definition
        isInTbl = verifyEpoch(tbl, epoch);
        tSt = tbl.start_time(strcmp(tbl.tags, epoch));
        tSp = tbl.stop_time(strcmp(tbl.tags, epoch));
        
        % Find timestamp data that falls within that interval and recover
        % their sample-indices (this is necessary because the data may not
        % be continuous, but rather have time-skips, which would be seen in
        % the timestamps
        tstmpAll = elecSeries.timestamps.load;
        samps(1) = find((tstmpAll >= tSt), 1, 'first');
        samps(2) = find((tstmpAll < tSp), 1, 'last');
        
        
    otherwise
        error('epoch string not recognized!')
        
end
        


%% Determine which channel of lfp to load
% if not specified, load all

% If channel not specified, ch indicates spread of all channels
if isempty(ch), ch = [1, nChans]; end



%% Final loading

lfp = elecSeries.data.load([samps(1), ch(1)], [samps(end), ch(end)]);
timestamps = elecSeries.timestamps.load(samps(1), samps(end));
fs = elecSeries.starting_time_rate; % Hz



end

%% SUB-FUNCTIONS

function TF = verifyEpoch(tbl, epochStr)
TF = false;

if ~any(strcmp(tbl.tags, epochStr))
    error('Specified epoch tag not found in nwb epoch intervals data.');
    
else
    TF = true;

end

end