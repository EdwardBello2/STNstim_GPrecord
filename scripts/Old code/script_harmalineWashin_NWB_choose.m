% batch for generating lots of iterations of the washin analysis

% NOTE: this script will only work on harmaline-day 

% TO-DO
% - test script as it currently is
% - test ability to Cache
% - refactor

%% initialize "cfg" struct with settings to be passed thru the workflow
% projConfig



%% SET variables

% get subselection from table for a given session
sessionMetaData.dose = '6 mg/kg'; % str, '2 mg/kg' | '4 mg/kg' | '6 mg/kg' | '8 mg/kg'
% sessionMetaData.exposureNum = '1st'; % str, '1st' | '2nd'
sessionMetaData.sessionType = 'exp'; % str, 'baseline' | 'exp'

% exp day epochs
% sessionMetaData.epochs = {'naivePreExp', 'naiveDbsExp', 'harWashinExp', ...
%     'harTaskExp', 'harDbsExp', 'harWashoutExp'};

sessionMetaData.epochs = {'naivePreExp', 'harWashinExp', 'harWashoutExp'};


LFP_BAND = [12, 35]; % Hz
CHPAIR = 2;
movingwin = [3, 1]; % seconds, same parameter as for mtspecgramc



%% High-level function for getting lfp power over time, normalized
clear bnPow times
% 1st exposure
sessionMetaData.exposureNum = '1st'; % str, '1st' | '2nd'
nwb = nwbRead_specSess(cfg, sessionMetaData); % get specific nwb file

[bnPow{1}, times{1}] = func_harmalineWashin_NWB_spec(...
    nwb, movingwin, cfg, LFP_BAND, CHPAIR, sessionMetaData);


% 2nd exposure
sessionMetaData.exposureNum = '2nd'; % str, '1st' | '2nd'
nwb = nwbRead_specSess(cfg, sessionMetaData); % get specific nwb file

[bnPow{2}, times{2}] = func_harmalineWashin_NWB_spec(...
    nwb, movingwin, cfg, LFP_BAND, CHPAIR, sessionMetaData);





%% PLOT the results of function

figure; ax = axes;
plot(times{1}, bnPow{1}, times{2}, bnPow{2});
% ylabel('Norm to naivePre, A.U.');
ylabel('10log Power');
xlabel('Time (min)');
legend('1st', '2nd')


% Craft title:
DOSE = sessionMetaData.dose;
EXPOS = sessionMetaData.exposureNum;
SESSTYPE = sessionMetaData.sessionType;
for i = 1:8, cStr(i) = {['C' num2str(i-1)]}; end
for i = 1:7, cpStr(i) = {[cStr{i} '-' cStr{i+1}]}; end

titStr = [cpStr{CHPAIR} ' ' DOSE ' peri-harmaline [' num2str(LFP_BAND(1)) ' - ' ...
    num2str(LFP_BAND(2)) ' Hz] band-power'];
title(titStr);

% titStr = [cpStr{CHPAIR} ' ' EXPOS ' ' DOSE ' peri-harmaline [' num2str(LFP_BAND(1)) ' - ' ...
%     num2str(LFP_BAND(2)) ' Hz] band-power'];
% title(titStr);


