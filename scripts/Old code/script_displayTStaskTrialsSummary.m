% script for reading in summary data of touchscreen task for a given session
% Input is the text content of the txt files output by
% BasicParadigmEnhanced custom software from the Vitek group that I used
% for the harmaline study in Uva


% File to be read in:
fn = 'MonkeyX9-23-2019-11-50-15-AM_2';

pn = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TStask BasicParadigmEnhanced\harmalineDoseStudy\';

inChar = fileread([pn fn '\' fn '.txt']);

[Summary, ReachTrials] = parseInfoRound(inChar);

open ReachTrials

%%
trialNumber = str2double(regexp(inChar, 'TrialNumbe:\s+\d+', 'match'));
trialNumber = str2double(regexp(inChar, 'TrialNumbe:\s+(\d+)', 'match'));

TriggersPattern = 'TrialNumber:\s+(\d+)';
data = regexp(inChar, TriggersPattern, 'tokens');
nTokens = numel(data);
unwrapd = zeros(nTokens, 1);
for iToken = 1:numel(data)
    unwrapd(iToken,1) = str2double(data{1,iToken}{1,1});
    
end
trialNumber = unwrapd;