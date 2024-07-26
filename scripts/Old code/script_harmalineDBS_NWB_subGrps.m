% script for midway point beween choosing just one combo and batch-saving
% images of every possible combo

% NOTE: this script will only work on harmaline-day 


% THIS SCRIPT WILL PLOT ALL 8 experiments for the parameters described
% below.

% TO-DO


%% initialize "cfg" struct with settings to be passed thru the workflow
% projConfig

figPos{1} = [2344 560 560 420];
figPos{2} = [2905 560 560 420];
figPos{3} = [2343 55 560 420];
figPos{4} = [2905 55 560 420];


%% SET variables

CHPAIR = [5, 6];


% get subselection from table for a given session
DOSES = {'2 mg/kg', '4 mg/kg', '6 mg/kg', '8 mg/kg'};
% sessionMetaData.dose = DOSES
EXPOS = {'1st', '2nd'};

EPOCHS = {'naiveDbsExp', 'harDbsExp'};

% sessionMetaData.exposureNum = '1st'; % str, '1st' | '2nd'
sessionMetaData.sessionType = 'exp'; % str, 'baseline' | 'exp'

% exp day epochs
% sessionMetaData.epochs = {'naivePreExp', 'naiveDbsExp', 'harWashinExp', ...
%     'harTaskExp', 'harDbsExp', 'harWashoutExp'};

% sessionMetaData.epochs = {'naiveDbsExp'};


% LFP_BAND = [2, 4]; % Hz
movingwin = [3, 1]; % seconds, same parameter as for mtspecgramc

% Run analysis and plotting for each experiment
fg = figure; 
fg.Position = [1940 137 1841 835];



% Test case of just one nwb file:


% sessionMetaData.dose = '2 mg/kg';
% sessionMetaData.exposureNum = '1st';
% sessionMetaData.sessionType = 'exp'; % str, 'baseline' | 'exp'
% 
% 
% nwb = nwbRead_specSess(cfg, sessionMetaData); % get specific nwb file
% 
% [S,tspec,f,tst,~] = specgram_nwblfp(nwb, movingwin, cfg, ...
%     'epoch', sessionMetaData.epochs{1}, ...
%     'channel', CHPAIR);
% 
% ax1 = subplot(2, 4, 1);
% ax1.Parent = fg;
% hold(ax1, 'on');
% surf(tspec, f, pow2db(S.^2)', ...
%     'Parent', ax1, ...
%     'EdgeColor', 'none');
% 
% % Create ylabel
% ylabel('Frequency (Hz)');
% 
% % Create xlabel
% xlabel('Time (s)');
% 
% grid(ax1,'on');
% axis(ax1,'tight');
% colorbar('peer',ax1);
% 
% 
% 
% 
% sessionMetaData.dose = '2 mg/kg';
% sessionMetaData.exposureNum = '2nd';
% sessionMetaData.sessionType = 'exp'; % str, 'baseline' | 'exp'
% 
% 
% nwb = nwbRead_specSess(cfg, sessionMetaData); % get specific nwb file
% 
% [S,tspec,f,tst,~] = specgram_nwblfp(nwb, movingwin, cfg, ...
%     'epoch', sessionMetaData.epochs{1}, ...
%     'channel', CHPAIR);
% 
% ax2 = subplot(2, 4, 5);
% ax2.Parent = fg;
% hold(ax2, 'on');
% surf(tspec, f, pow2db(S.^2)', ...
%     'Parent', ax2, ...
%     'EdgeColor', 'none');
% 
% % Create ylabel
% ylabel('Frequency (Hz)');
% 
% % Create xlabel
% xlabel('Time (s)');
% 
% grid(ax1,'on');
% axis(ax1,'tight');
% colorbar('peer',ax2);
% 
% 
% 
% 
% dispSpecgram(tspec, f, S')
% dispSpecgram(tspec, f, (S.^2)')
% dispSpecgram(tspec, f, pow2db(S.^2)')


nDoses = numel(DOSES);
nExp = numel(EXPOS);
nEpochs = numel(EPOCHS);
tic

i = 1; % index for subplot number
for iD = 1:nDoses
    for iEx = 1:nExp
        for iEp = 1:nEpochs
            sessionMetaData.dose = DOSES{iD};
            sessionMetaData.exposureNum = EXPOS{iEx};
            sessionMetaData.epochs = EPOCHS{iEp};
            sessionMetaData.sessionType = 'exp'; % str, 'baseline' | 'exp'


            nwb = nwbRead_specSess(cfg, sessionMetaData); % get specific nwb file

            [S,tspec,f,tst,~] = specgram_nwblfp(nwb, movingwin, cfg, ...
                'epoch', sessionMetaData.epochs, ...
                'channel', CHPAIR);

            ax(i) = subplot(2, 8, i);
            ax(i).Parent = fg;
            hold(ax(i), 'on');
            surf(tspec, f, pow2db(S.^2)', ...
                'Parent', ax(i), ...
                'EdgeColor', 'none');

            % Create ylabel
%             ylabel('Frequency (Hz)');

            % Create xlabel
            xlabel('Time (s)');

            grid(ax(i),'on');
            axis(ax(i),'tight');
%             colorbar(ax(i));  

            titStr = [DOSES{iD} ', ' EXPOS{iEx} ', ' EPOCHS{iEp}];
            title(titStr);

            i = i+1; % increment to next axes position in subplot
            
        end
        
    end
    
end
toc
linkaxes(ax)
set(ax(:), 'YLim', [0, 150]);
set(ax(:),'CLim',[-280 -220]);

% 
% for iD = 1:nDoses
%     % High-level function for getting lfp power over time, normalized
%     clear bnPow times
%     
%     % Set harmaline dose to plot for this iteration
%     sessionMetaData.dose = DOSES{iD};
%     
%     
%     % 1st exposure
%     sessionMetaData.exposureNum = '1st'; % str, '1st' | '2nd'
% %     tic
%     nwb = nwbRead_specSess(cfg, sessionMetaData); % get specific nwb file
% %     toc
% 
% %     tic
%     [bnPow{1}, times{1}] = func_harmalineWashin_NWB_spec(...
%         nwb, movingwin, cfg, LFP_BAND, CHPAIR, sessionMetaData);
% %     toc
% 
%     % 2nd exposure
%     sessionMetaData.exposureNum = '2nd'; % str, '1st' | '2nd'
% %     tic
%     nwb = nwbRead_specSess(cfg, sessionMetaData); % get specific nwb file
% %     toc
%     
% %     tic
%     [bnPow{2}, times{2}] = func_harmalineWashin_NWB_spec(...
%         nwb, movingwin, cfg, LFP_BAND, CHPAIR, sessionMetaData);
% %     toc
%     
% %     tic
%     % PLOT the results of function
%     ax(iD) = subplot(2,2,iD);
%     
%     plot(times{1}, bnPow{1}, times{2}, bnPow{2});
%     % ylabel('Norm to naivePre, A.U.');
%     ylabel('10log Power');
%     xlabel('Time (min)');
%     legend('1st', '2nd')
% 
% 
%     % Craft title:
%     DOSE = sessionMetaData.dose;
%     EXPOS = sessionMetaData.exposureNum;
%     SESSTYPE = sessionMetaData.sessionType;
%     for i = 1:8, cStr(i) = {['C' num2str(i-1)]}; end
%     for i = 1:7, cpStr(i) = {[cStr{i} '-' cStr{i+1}]}; end
% 
%     titStr = [cpStr{CHPAIR} ' ' DOSE ' peri-harmaline [' num2str(LFP_BAND(1)) ' - ' ...
%         num2str(LFP_BAND(2)) ' Hz] band-power'];
%     title(titStr);
%     toc
%     
%     
%     
%     
%     
% end
% linkaxes(ax); % allows me to easily manually change the zoom on all at once
% 
% toc
% 
% hallelujah()
% 

