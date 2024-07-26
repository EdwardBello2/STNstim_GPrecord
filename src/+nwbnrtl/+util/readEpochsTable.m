function epoch_tbl = readEpochsTable(timeIntervals)

% element id
id = timeIntervals.id.data.load;

% fill epoch tags column
tag = timeIntervals.tags.data.load;

% fill start time column
start_time = timeIntervals.start_time.data.load;

% fill stop time column
stop_time = timeIntervals.stop_time.data.load;

epoch_tbl = [table(id), table(tag), table(start_time), table(stop_time)];

epoch_tbl = sortrows(epoch_tbl, 1);


end