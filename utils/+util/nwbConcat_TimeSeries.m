function tSeriesConc = nwbConcat_TimeSeries(nwbTSloadpath, nwbBlk)
% Create a concatenated TimeSeries object from multiple nwb files in one
% session. Assumes that the same TimeSeries data is present in all nwb
% files we're concatenating from


str = nwbTSloadpath;
strFinal = join(["tSeriesConc = nwbBlk(1)." str ";"], "");
eval(strFinal);


% Gather TimeSeries object for each NWB block in array: tSeriesBlk
nBlocks = length(nwbBlk);
for iBlk = 1:nBlocks
    strFinal = join(["tSeriesBlk(" num2str(iBlk) ") = nwbBlk(" num2str(iBlk) ")." str ";"], "");
    eval(strFinal);
    
end



%% Iteratively append data and timestamps

% Only things to update are "data" and "timestamps"
refDT = datetime(nwbBlk(1).session_start_time);

% FOR each nwb block:
data = [];
timeStamps = [];
for iBlk = 1:nBlocks

    dataBlk = tSeriesBlk(iBlk).data.load;
    
    % Get the re-referenced starting time
    refSec = seconds(datetime(nwbBlk(iBlk).timestamps_reference_time) - refDT);
    
    % Get (or generate) re-referenced timestamps
    startTime = tSeriesBlk(iBlk).starting_time;
    fs = tSeriesBlk(iBlk).starting_time_rate;
    tst = ((0:(length(dataBlk)-1)) * (1/fs)) + startTime;
    tstReref = tst + refSec;
    
    % Append timestamps and data for this iteration
    timeStamps = [timeStamps; tstReref'];
    data = [data; dataBlk];

end 



%% Save concatenated data into NWB session object

tSeriesConc.data = data;
tSeriesConc.timestamps = timeStamps;



end