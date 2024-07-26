% script for reading one row of data from a .txt generated by
% BasicParadigmEnhanced


% File to be read in:
fn = 'MonkeyX4-7-2019-9-31-46-AM_3';


% Read in text file to one long character vector
pn = 'L:\Shared drives\Johnson\TouchscreenTask\TrainingRecords\Uva';
inChar = fileread([pn, '\', fn, '\', fn, '.txt.']);

% Convert char vector into table row of data
Row = parseInfoRound(inChar);


% Append info to the row
filename = {fn};
Fn = table(filename);

Round = 1; % default round assignment for any given round is 1; user must change manually

Rnd = table(Round);

RowFinal = [Row(:,1:3), Rnd, Fn, Row(:,4:end)];

% writetable(RowFinal)
% writetable(RowFinal, 'test', 'FileType', 'spreadsheet')

