function tSeriesConc = concat_TimeSeries(tSeriesBlk, session_start_time)
% Create a concatenated TimeSeries object from multiple nwb files in one
% session. Assumes that the same TimeSeries data is present in all nwb
% files we're concatenating from


% str = nwbTSloadpath;
% strFinal = join(["tSeriesConc = nwbBlk(1)." str ";"], "");
% eval(strFinal);
if length(tSeriesBlk) ~= length(session_start_time)
    error('Both inputs need to be same length!')
    
end

% Initialize the concatenated TimeSeries based on first original one
tSeriesConc = tSeriesBlk(1);

% Gather TimeSeries object for each NWB block in array: tSeriesBlk
nBlocks = length(tSeriesBlk);
% for iBlk = 1:nBlocks
%     strFinal = join(["tSeriesBlk(" num2str(iBlk) ") = nwbBlk(" num2str(iBlk) ")." str ";"], "");
%     eval(strFinal);
%     
% end



%% Iteratively append data and timestamps

% Only things to update are "data" and "timestamps"
refDT = datetime(session_start_time(1)); % based on first datetime object

% FOR each nwb block:
data = [];
timeStamps = [];
for iBlk = 1:nBlocks

    dataBlk = tSeriesBlk(iBlk).data.load;
    
    % Get the re-referenced starting time
    refSec = seconds(datetime(session_start_time(iBlk)) - refDT);
    
    % Get (or generate) re-referenced timestamps
    switch tSeriesConc.data_continuity
        case {'instantaneous', 'step'}
            tst = tSeriesBlk(iBlk).timestamps.load;
            
        case 'continuous'
            startTime = tSeriesBlk(iBlk).starting_time;
            fs = tSeriesBlk(iBlk).starting_time_rate;
            tst = ((0:(length(dataBlk)-1)) * (1/fs)) + startTime;
            tst = tst';
            
        otherwise
            error('bla!!!')
            
    end
    tstReref = tst + refSec;

    
    % Append timestamps and data for this iteration
    timeStamps = [timeStamps; tstReref];
    data = [data; dataBlk];

end 



%% Save concatenated data into NWB session object

tSeriesConc.data = data;
tSeriesConc.timestamps = timeStamps;



end