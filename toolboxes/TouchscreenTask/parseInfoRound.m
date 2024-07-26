function [EndSummary, ReachTrials] = parseInfoRound(inText)
% Reads in a long char vector, specificially the output .txt generated by
% BasicParadigmEnhanced program during a touchscreen task recording.
% Creates a table-row of data about the entire recording. 

% Author: Ed Bello
% Created: 2019/03/29
% To-Do:
% - Find a way to account for those files where number of data for various
% fields seems to differ. Example: 85 reaction times detected and 87 reach
% times detected. No idea why this happens. For now output things in
% ReachTrials struct rather than a table, since tables require all columns
% to have same number of elements. Strictly speaking, the struct doesn't
% for sure match data points across its fields with same trials data.



TriggersPattern = 'Triggers:\s+(\d+)';
data = regexp(inText, TriggersPattern, 'tokens');
Triggers = str2double(cell2mat(data{1}));

PresentationsPattern = 'Presentations:\s+(\d+)';
data = regexp(inText, PresentationsPattern, 'tokens');
Presentations = str2double(cell2mat(data{1}));

SuccessTrialsPattern = 'SuccessTrials:\s+(\d+)';
data = regexp(inText, SuccessTrialsPattern, 'tokens');
SuccessTrials = str2double(cell2mat(data{1}));

ErrorTrialsPattern = 'ErrorTrials:\s+(\d+)';
data = regexp(inText, ErrorTrialsPattern, 'tokens');
ErrorTrials = str2double(cell2mat(data{1}));

AbortedTrialsPattern = 'AbortedTrials:\s+(\d+)';
data = regexp(inText, AbortedTrialsPattern, 'tokens');
AbortedTrials = str2double(cell2mat(data{1}));

TargetHoldFailureTrialsPattern = 'TargetHoldFailureTrials:\s+(\d+)';
data = regexp(inText, TargetHoldFailureTrialsPattern, 'tokens');
TargetHoldFailureTrials = str2double(cell2mat(data{1}));

ReplicationsCompletedPattern = 'ReplicationsCompleted:\s+(\d+)';
data = regexp(inText, ReplicationsCompletedPattern, 'tokens');
ReplicationsCompleted = str2double(cell2mat(data{1}));

% Get total time of this round
tBegPattern = 'Date:\s+\d+/\d+/\d+\s+(\d+:\d+:\d+)';

monthPattern = 'Date:\s+(\d+)';
data = regexp(inText, monthPattern, 'tokens');
MM = str2double(cell2mat(data{1}));

dayPattern = 'Date:\s+\d+/(\d+)';
data = regexp(inText, dayPattern, 'tokens');
DD = str2double(cell2mat(data{1}));

yearPattern = 'Date:\s+\d+/\d+/(\d+)';
data = regexp(inText, yearPattern, 'tokens');
YYYY = str2double(cell2mat(data{1}));



EndSummary = table(YYYY, MM, DD, ...
          Triggers, ...
          Presentations, ...
          SuccessTrials, ...
          ErrorTrials, ...
          AbortedTrials, ...
          TargetHoldFailureTrials, ...
          ReplicationsCompleted);
      
      
%% Gathering data related to each individual reach trial

       ReachTrials.TrialNumber = unwrapTokensNumerical(regexp(inText, 'TrialNumber:\s+(\d+)', 'tokens'));
    ReachTrials.TriggerNumber = unwrapTokensNumerical(regexp(inText, 'TriggerNumber:\s+(\d+)', 'tokens'));
  ReachTrials.BaselineDuration = unwrapTokensNumerical(regexp(inText, 'BaselineDuration:\s+(\d+.\d+)', 'tokens'));
         ReachTrials.SizeIndex = unwrapTokensNumerical(regexp(inText, 'SizeIndex:\s+(\d+)', 'tokens'));
     ReachTrials.PositionIndex = unwrapTokensNumerical(regexp(inText, 'PositionIndex:\s+(\d+)', 'tokens'));
ReachTrials.InstructionalDelay = unwrapTokensNumerical(regexp(inText, 'InstructionalDelay:\s+(\d+)', 'tokens'));
      ReachTrials.GoCueWrtBase = unwrapTokensNumerical(regexp(inText, 'GoCueWrtBase:\s+(\d+.\d+)', 'tokens'));
      ReachTrials.ReactionTime = unwrapTokensNumerical(regexp(inText, 'ReactionTime:\s+(\d+.\d+)', 'tokens'));
         ReachTrials.ReachTime = unwrapTokensNumerical(regexp(inText, 'ReachTime:\s+(\d+.\d+)', 'tokens'));
%  CoordinateTouched_x = regexp(inText, 'CoordinateTouched:\s+(\d+),\d+', 'tokens')
%  CoordinateTouched_y = regexp(inText, 'CoordinateTouched:\s+\d+,(\d+)', 'tokens')
        ReachTrials.TrialSuccess = unwrapTokensTF(regexp(inText, 'TrialSuccess:\s+(\w+)', 'tokens'));


%         
% ReachTrials = table(TrialNumber, TriggerNumber, BaselineDuration, ...
%     SizeIndex, PositionIndex, InstructionalDelay, GoCueWrtBase, ...
%     ReactionTime, ReachTime, TrialSuccess);



end

%% SUB-FUNCTIONS

function unwrapd = unwrapTokensNumerical(data)
nTokens = numel(data);
unwrapd = zeros(nTokens, 1);
for iToken = 1:numel(data)
    unwrapd(iToken,1) = str2double(data{1,iToken}{1,1});
    
end

end

function unwrapd = unwrapTokensTF(data)
nTokens = numel(data);
unwrapd = false(nTokens, 1);
for iToken = 1:numel(data)
    tfStr = (data{1,iToken}{1,1});
    if strcmp(tfStr, 'True')
        unwrapd(iToken,1) = true; 
        
    elseif strcmp(tfStr, 'False')
        unwrapd(iToken,1) = false;
        
    else
        error('Unexpected string for TF string')
        
    end
    
end

end

