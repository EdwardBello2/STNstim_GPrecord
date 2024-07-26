% script for reading all rows of each round, compiled into one big xlsx
% sheet.


% % File to be read in:
% fn = 'MonkeyX3-1-2019-11-46-54-AM_5';

pn = 'L:\Shared drives\Johnson\TouchscreenTask\TrainingRecords\Uva';

% Get a cell array of all rounds folders recorded so far
contents = struct2cell(dir(pn));
names = contents(1,:);
names(1:2) = [];
names(end) = [];

nRounds = numel(names);

TS = array2table(zeros(0,12));
TS.Properties.VariableNames = {'YYYY', 'MM', 'DD', 'Round', 'Filename', ...
    'Triggers','Presentations', ...
    'SuccessTrials', 'ErrorTrials', 'AbortedTrials', ...
    'TargetHoldFailureTrials', 'ReplicationsCompleted'};

for iRound = 1:nRounds
    fn = names{iRound};

    % Read in text file to one long character vector
    inChar = fileread([pn, '\', fn, '\', fn, '.txt']);

    % Convert char vector into table row of data
    try
        Row = parseInfoRound(inChar);

        % Append info to the row
        Filename = {fn};
        Fn = table(Filename);
        Round = 1;
        Rnd = table(Round);

        RowFinal = [Row(:,1:3), Rnd, Fn, Row(:,4:end)];

        % add the row to the overall table:
        TS = [TS; RowFinal];
        
    catch
        disp(['text file could not be parsed: ', fn]);
        
    end

end

%%
% writetable(TS, 'test', 'FileType', 'spreadsheet')




