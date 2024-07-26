function timestamps = readTimeSeriesTimestamps(timeseries, varargin)
%   D = LOADTIMESERIESDATA(TIMESERIES) loads all the data from the timestamps field
%   of the timeseries object. D is of shape samples x 1.
%
%   D = LOADTIMESERIESDATA(TIMESERIES, INTERVALTIMES) loads timestamps field
%   of the timeseries object, between the range specified by intervalTimes
%   (seconds). Times are relative to session start time synchronization.
%   Defaults to using the timestamps values if they exist, otherwise just
%   uses start time and sampling frequency. 
%
%   D = LOADTIMESERIESDATA(TIMESERIES, INTERVALTIMES, DOWNSAMPLE_FACTOR)
%   specifies a temporal downsampling for D. Default is 1.
%   
%   D = LOADTIMESERIESDATA(TIMESERIES, ___ , Name, Value)
%
%   'intervalSamples' -- Sample-Range
%   [] (default) | two-element array of positive integers
%   specifying range of samples to
%   subselect from data, using sample instead of time in seconds. Values
%   specified in this property will override whatever input was specified
%   for INTERVALTIMES.
%   
%


%% Input Parser

% Load only data within a certain epoch
p = inputParser;

addRequired(p, 'timeseries');
addOptional(p, 'intervalTimes', [0, Inf]);
addOptional(p, 'downsample_factor', 1);
% addOptional(p, 'electrode', []);
addParameter(p, 'intervalSamples', []);
parse(p, timeseries, varargin{:});


intervalTimes = p.Results.intervalTimes;
downsample_factor = p.Results.downsample_factor;
% electrode = p.Results.electrode;
intervalSamples = p.Results.intervalSamples;



%%

if isa(timeseries.data, 'types.untyped.DataPipe')
    dims = timeseries.data.internal.dims;
        
else % usual case for normal TimeSeries objects
    dims = timeseries.data.dims;

end


if ~exist('intervalTimes','var') || isempty(intervalTimes)
    intervalTimes = [0 Inf];
    
end

% STRIDE = ones(1, length(dims));
if ~exist('downsample_factor','var') || isempty(downsample_factor)
    STRIDE(1) = 1;
    
else
    STRIDE(1) = downsample_factor;
    
end

% if ~exist('electrode', 'var')
%     electrode = [];
%     
% end



%%
    
% else % use sampling rate
    if intervalTimes(1)
        if ~isempty(timeseries.timestamps)
            disp('Searching timestamps for interval...');
            tst = timeseries.timestamps.load;
            disp('Done Searching');
            start_ind = fastsearch(tst, intervalTimes(1), 1);

        else % use sampling rate
            fs = timeseries.starting_time_rate;
            t0 = timeseries.starting_time;
            if intervalTimes(1) < t0
                error('interval bounds outside of time range');

            end
            start_ind = (intervalTimes(1) - t0) * fs;

        end

    else
        start_ind = 1;

    end

    if isfinite(intervalTimes(2))        
        if ~isempty(timeseries.timestamps)
            end_ind = fastsearch(tst, intervalTimes(2), -1);

        else
            fs = timeseries.starting_time_rate;
            t0 = timeseries.starting_time;
            if intervalTimes(2) > (dims(end) * fs + t0)
                error('interval bounds outside of time range');

            end
            end_ind = (intervalTimes(2) - t0) * fs;

        end

    else
        end_ind = dims(1);

    end

% end



%% Load Data depending on which electrode is specified

% START = ones(1, length(dims));
START = start_ind;

% END = ones(1, length(dims));
END = end_ind;

nSamps = end_ind - start_ind + 1;


% Special overriding case if the parameter "intervalSamples" was specified
if ~isempty(intervalSamples)
    START(1) = intervalSamples(1);
    END(1) = intervalSamples(2);
    nSamps = intervalSamples(2) - intervalSamples(1) + 1;

end


% % fill the data field based on the electrode input
% % disp('Loading timeseries data...')
% if length(electrode) > 1 % then we call the code recursively for each
%     data = NaN(nSamps, length(electrode));
%     for i = 1:length(electrode)
%         data(:,i) = nwbnrtl.util.loadTimeSeriesData(timeseries, [], ...
%             downsample_factor, electrode(i), ...
%             'intervalSamples', [START(1), END(1)]);
%         disp(['... electrode ' num2str(electrode(i)) ' COMPLETE'])
% 
%         
%     end
%     
% else
%     if isempty(electrode)
%         START(end) = 1;
%         END(end) = dims(end);
% 
%     else
%         START(end) = electrode;
%         END(end) = electrode;
% 
%     end
%     
%     disp('Loading electrode...')
    if isa(timeseries.timestamps, 'types.untyped.DataPipe')
        timestamps = timeseries.timestamps.internal.load(START, STRIDE, END);

    else % usual case for normal TimeSeries objects
        timestamps = timeseries.timestamps.load(START, STRIDE, END);

    end   
%     disp('Done loading all !')

    % Scale data properly
%     data = data * timeseries.data_conversion;

% end