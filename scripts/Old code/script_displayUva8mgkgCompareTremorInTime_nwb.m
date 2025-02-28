% script for comparing tremor power and timecourse between Mela, Kiwi, and
% Uva. 
clear; clc

%%
% Loading and saving pathways
clear nwb
nwbDataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';
projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';

% % Uva dose study data of choice:
% nwbFullPath = nwbDataPath;
% % SESSION_ID = 'TremoLfpDBS-190927-100155'; % 2 mg/kg day
% % SESSION_ID = 'TremoLfpDBS-191004-100637'; % 4 mg/kg day
% % SESSION_ID = 'TremoLfpDBS-191011-104322'; % 6 mg/kg day
% % SESSION_ID = 'TremoLfpDBS-191018-100615'; % 8 mg/kg day
% % SESSION_ID = 'TremoLfpDBS-191025-104651'; % 2 mg/kg day
% % SESSION_ID = 'TremoLfpDBS-191101-101430'; % 4 mg/kg day
% % SESSION_ID = 'TremoLfpDBS-191108-101829'; % 6 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191115-100127'; % 8 mg/kg day

% Remove Chronux from Matlab's searchpath to prevent a problem with
% jacknife function used by bootci
rmpath(genpath('C:\Users\bello043\Documents\GitHub\Thalamic-DBS-HarmalineTremor\toolboxes\chronux_2_10'))



%%

% CONFIGURATION & CONSTANTS
% decide whether to show all figures generated by workflow
cfg.displayWorkflowFigs = false;

% (the the rest of the constants See inside "tremorWorkflow" function at bottom of script)

cfg.fcDetrend = 1; % Hz

% Spectrogram parameters:
cfg.TWIN = 1; % seconds, sliding time-window of spectrogram
cfg.NOVERLAP = 0; % samples, no overlaps for these spectrogam windows

% Movement artifact removal
cfg.MADthresh = 3; % median limit for outliers

% Tremor frequency band
% trPeak = 10; % Hz, major tremor component
% trRange = 3; % limits of tremor band around major component
% cfg.tremorBand = [6, 12]; % Hz

% General movement frequency band
generalBand = [17, 21];

% Tremor power summary statistic & bootstrapped estimation of confidence interval
% myStatFun = @(x)mean(x);
cfg.NBOOTS = 1000;

% If any time-series recordings MUCH larger than the typical small size,
% break into smaller segments for display
cfg.tooBigTime = 240; % seconds
cfg.subSegTime = 90; % seconds
cfg.tooSmallTime = 30; % seconds

cfg.tooSparseData = 30;

% Baseline estimate statistic

% Test for segments vs baseline
statTest = @(x)signrank(x);

% False Discovery Rate for Benjamini-Hochberg procedure for multiple
% comparisons
fdr = 0.005;




% %% Mela's data
% 
% % Mela Pilot data of choice: 12 mg/kg day
% nwbFullPath = [nwbDataPath 'MelaPilot\'];
% % SESSION_ID = 'harPilotMela-150928-094900'; % 2 mg/kg day
% % SESSION_ID = 'harPilotMela-150930-082000'; % 6 mg/kg day
% SESSION_ID = 'harPilotMela-151013-072100'; % 12 mg/kg day
% 
% % Load nwb master file for this session, get acc TimeSeries object
% disp('Reading NWB file...')
% nwbM = nwbRead([nwbFullPath SESSION_ID '_MASTER.nwb']);
% disp('DONE!')
% 
% cfg.tremorBand = [8, 12]; % Hz
% % Calculate spectrogram-based tremor power for segments of recorded data
% [tremPowM, tstSegHarRefM] = func_gatherACCproc_fbandTimeSegs_har(nwbM, cfg);
% 
% 
% % Prepare that data into format that can be plotted by errorbar function
% [tmidM,powM,cinegM,ciposM,tnegM,tposM] = func_prepACCproc_fbandTimeSegs_errorbar(...
%     tremPowM, tstSegHarRefM, cfg);
% 
% 
% 
% %% Kiwi's data
% 
% % Kiwi Pilot data of choice: 10 mg/kg day
% nwbFullPath = [nwbDataPath 'KiwiPilot\'];
% SESSION_ID = 'harPilotKiwi-160517-122121';
% 
% % Load nwb master file for this session, get acc TimeSeries object
% disp('Reading NWB file...')
% nwbK = nwbRead([nwbFullPath SESSION_ID '_MASTER.nwb']);
% disp('DONE!')
% 
% cfg.tremorBand = [3, 8]; % Hz
% % Calculate spectrogram-based tremor power for segments of recorded data
% [tremPowK, tstSegHarRefK] = func_gatherACCproc_fbandTimeSegs_har(nwbK, cfg);
% 
% 
% % Prepare that data into format that can be plotted by errorbar function
% [tmidK,powK,cinegK,ciposK,tnegK,tposK] = func_prepACCproc_fbandTimeSegs_errorbar(...
%     tremPowK, tstSegHarRefK, cfg);
% 


%% Uva's first 8 mg/kg dose

% Uva Pilot data of choice: 8 mg/kg day
nwbFullPath = [nwbDataPath 'UvaPilot\'];
% SESSION_ID = 'harPilotUva-180219-132618'; % 8 mg/kg day
SESSION_ID = 'harPilotUva-180221-113120'; % 12 mg/kg day
% SESSION_ID = 'harPilotUva-180305-103016'; % 12 mg/kg day

% Load nwb master file for this session, get acc TimeSeries object
disp('Reading NWB file...')
nwbU1 = nwbRead([nwbFullPath SESSION_ID '_MASTER.nwb']);
disp('DONE!')

cfg.tremorBand = [8, 12]; % Hz
% Calculate spectrogram-based tremor power for segments of recorded data
[tremPowU1, tstSegHarRefU1] = func_gatherACCproc_fbandTimeSegs_har(nwbU1, cfg);

% Calculate spectrogram-based notremor power for segments of recorded data
cfg.tremorBand = generalBand; % Hz
[notremPowU1, ~] = func_gatherACCproc_fbandTimeSegs_har(nwbU1, cfg);

[powCompU1] = comparePow(tremPowU1, notremPowU1);




% Prepare that data into format that can be plotted by errorbar function
[tmidU1,powU1,cinegU1,ciposU1,tnegU1,tposU1] = func_prepACCproc_fbandTimeSegs_errorbar(...
    tremPowU1, tstSegHarRefU1, cfg);
[tmidU1,powU1,cinegU1,ciposU1,tnegU1,tposU1] = func_prepACCproc_fbandTimeSegs_errorbar(...
    powCompU1, tstSegHarRefU1, cfg);



%% Uva's second 

% Uva dose study data of choice:
nwbFullPath = nwbDataPath;
% SESSION_ID = 'TremoLfpDBS-190927-100155'; % 2 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191004-100637'; % 4 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191011-104322'; % 6 mg/kg day
SESSION_ID = 'TremoLfpDBS-191018-100615'; % 8 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191025-104651'; % 2 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191101-101430'; % 4 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191108-101829'; % 6 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191115-100127'; % 8 mg/kg day

% Load nwb master file for this session, get acc TimeSeries object
disp('Reading NWB file...')
nwbU2 = nwbRead([nwbFullPath SESSION_ID '_MASTER.nwb']);
disp('DONE!')

cfg.tremorBand = [8, 12]; % Hz
% Calculate spectrogram-based tremor power for segments of recorded data
[tremPowU2, tstSegHarRefU2] = func_gatherACCproc_fbandTimeSegs_har(nwbU2, cfg);

% Calculate spectrogram-based notremor power for segments of recorded data
cfg.tremorBand = generalBand; % Hz
[notremPowU2, ~] = func_gatherACCproc_fbandTimeSegs_har(nwbU2, cfg);

[powCompU2] = comparePow(tremPowU2, notremPowU2);


% Prepare that data into format that can be plotted by errorbar function
[tmidU2,powU2,cinegU2,ciposU2,tnegU2,tposU2] = func_prepACCproc_fbandTimeSegs_errorbar(...
    tremPowU2, tstSegHarRefU2, cfg);
[tmidU2,powU2,cinegU2,ciposU2,tnegU2,tposU2] = func_prepACCproc_fbandTimeSegs_errorbar(...
    powCompU2, tstSegHarRefU2, cfg);



%% Uva's third

% Uva dose study data of choice:
nwbFullPath = nwbDataPath;
% SESSION_ID = 'TremoLfpDBS-190927-100155'; % 2 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191004-100637'; % 4 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191011-104322'; % 6 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191018-100615'; % 8 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191025-104651'; % 2 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191101-101430'; % 4 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191108-101829'; % 6 mg/kg day
SESSION_ID = 'TremoLfpDBS-191115-100127'; % 8 mg/kg day

% Load nwb master file for this session, get acc TimeSeries object
disp('Reading NWB file...')
nwbU3 = nwbRead([nwbFullPath SESSION_ID '_MASTER.nwb']);
disp('DONE!')

cfg.tremorBand = [8, 12]; % Hz
% Calculate spectrogram-based tremor power for segments of recorded data
[tremPowU3, tstSegHarRefU3] = func_gatherACCproc_fbandTimeSegs_har(nwbU3, cfg);

% Calculate spectrogram-based notremor power for segments of recorded data
cfg.tremorBand = generalBand; % Hz
[notremPowU3, ~] = func_gatherACCproc_fbandTimeSegs_har(nwbU3, cfg);

[powCompU3] = comparePow(tremPowU3, notremPowU3);


% Prepare that data into format that can be plotted by errorbar function
[tmidU3,powU3,cinegU3,ciposU3,tnegU3,tposU3] = func_prepACCproc_fbandTimeSegs_errorbar(...
    tremPowU3, tstSegHarRefU3, cfg);
[tmidU3,powU3,cinegU3,ciposU3,tnegU3,tposU3] = func_prepACCproc_fbandTimeSegs_errorbar(...
    powCompU3, tstSegHarRefU3, cfg);



%% Display All results

figure;
set(gcf, 'Position', [2083 66 1735 919]);

% First
s = 1;
ax(s) = subplot(3,1,s);
errorbar(tmidU1,powU1,cinegU1,ciposU1,tnegU1,tposU1, '.', ...
    'CapSize', 0, 'LineWidth', 3);
grid on; hold on;
ax(s).YLim(1) = 0;
plot([0 0], [ax(s).YLim], 'r')
ylabel('8-12Hz ACC Power (A.U.)')
legend('First', 'Location', 'northeastoutside')

title('Tremor Power following multiples Harmaline exposures to 8 mg/kg')


% Second
s = 2;
ax(s) = subplot(3,1,s);
errorbar(tmidU2,powU2,cinegU2,ciposU2,tnegU2,tposU2, '.', ...
    'CapSize', 0, 'LineWidth', 3);
grid on; hold on;
% ax(s).YLim = [0 0.25];
plot([0 0], [ax(s).YLim], 'r')
ylabel('8-12Hz ACC Power (A.U.)')
legend('Second', 'Location', 'northeastoutside')



% Third
s = 3;
ax(s) = subplot(3,1,s);
errorbar(tmidU3,powU3,cinegU3,ciposU3,tnegU3,tposU3, '.', ...
    'CapSize', 0, 'LineWidth', 3);
grid on; hold on;
% ax(s).YLim = [0 0.25];
plot([0 0], [ax(s).YLim], 'r')
ylabel('8-12Hz ACC Power (A.U.)')
legend('Third', 'Location', 'northeastoutside')

xlabel('Time since harmaline inject (mins)');
set(ax(:), 'XLim', [-40 280]);
set(ax(:), 'YLim', [0 3.2]);





%% SUB-FUNCTIONS

function [powComp] = comparePow(tremPowM, notremPowM)

nSegs = length(tremPowM);
for iSeg = 1:nSegs
    powComp{iSeg,1} = tremPowM{iSeg} ./ notremPowM{iSeg};
    
end

end

function [accSegFinal, tstSegFinal] = segmentacc(acc, tst, ctrl, maxSamp)
% separate acc data and timestamps into time segments based on recordings.
% Large recordings are further subdivided into smaller time segements according
% to max segment time. maxTime (seconds). ctrl has zero at the start of
% every recording. 

% for acc and tst composed of N separate recordings, ctrl tracks which data
% belongs to which N recording

% detect the number of segments
nSegs = max(ctrl);
for iSeg = 1:nSegs
    accSeg{iSeg,1} = acc((ctrl == iSeg),1:3);
    tstSeg{iSeg,1} = tst((ctrl == iSeg));
    nSamps(iSeg) = sum(ctrl == iSeg);
    
end

% If any segments are too long according to maxSamp 
accSegFinal = {};
tstSegFinal = {};
if exist('maxSamp', 'var')
    for iSeg = 1:nSegs
        if nSamps(iSeg) > maxSamp % Too big! subdivide the segment into smaller ones
            nSubs = floor(nSamps(iSeg)/maxSamp) + 1;
            tempAcc = accSeg{iSeg,1};
            tempTst = tstSeg{iSeg,1};
            
            % mark subdivisions with indices
            xb = 1:(nSubs-1);
            bIdx = [1, xb*maxSamp + 1];
            eIdx = [xb*maxSamp, nSamps(iSeg)];
            
            % loop-fill the new subsections to be added
            for iSub = 1:nSubs
                accSegFinal = [accSegFinal; {tempAcc(bIdx(iSub):eIdx(iSub),1:3)}];
                tstSegFinal = [tstSegFinal; {tempTst(bIdx(iSub):eIdx(iSub))}];
                
            end
            
%             xb = [1:nRowsAdd]';
% 
%             % begIdxCorr
%             begIdxSmall = begIdx(i) + xb*maxSamp + 1;
%             begIdxCorr = [begIdxCorr; begIdx(i); begIdxSmall];
% 
%             % endIdxCorr
%             endIdxSmall = begIdx(i) + xb*maxSamp;
%             endIdxCorr = [endIdxCorr; endIdxSmall; endIdx(i)];

        else
            accSegFinal = [accSegFinal; accSeg(iSeg,1)];
            tstSegFinal = [tstSegFinal; tstSeg(iSeg,1)];
        
        end
    
        
        
    end
    
else % simple case where maxSamp is not entered into the function, no subdivision done
    accSegFinal = accSeg;
    tstSegFinal = tstSeg;
    
end




% 
% 
% begIdx = find(ctrl == 0);
% endIdx = [(begIdx(2:end)-1); length(ctrl)];
% % nSamps = diff([begIdx endIdx], [], 2);
% 
% % enforce maximum samples for each segment
% % isToobig = nSamps > maxSamp;
% nSamps = diff([begIdx endIdx], [], 2);
% begIdxCorr = [];
% endIdxCorr = [];
% for i = 1:length(begIdx)
% 
%     if diff([begIdx(i) endIdx(i)], [], 2) > maxSamp % Too big! break up the segment into smaller ones
%         nRowsAdd = floor(nSamps(i)/maxSamp); 
%         xb = [1:nRowsAdd]';
%         
%         % begIdxCorr
%         begIdxSmall = begIdx(i) + xb*maxSamp + 1;
%         begIdxCorr = [begIdxCorr; begIdx(i); begIdxSmall];
%         
%         % endIdxCorr
%         endIdxSmall = begIdx(i) + xb*maxSamp;
%         endIdxCorr = [endIdxCorr; endIdxSmall; endIdx(i)];
%         
%     else
%         begIdxCorr = [begIdxCorr; begIdx(i)];
%         endIdxCorr = [endIdxCorr; endIdx(i)];
%         
%     end
%     
%     
%     
% end


% Gather data segments
% nSegs = length(begIdxCorr);
% for iSeg = 1:nSegs
%     accSeg{iSeg,1} = acc(begIdxCorr(iSeg):endIdxCorr(iSeg),:);
%     tstSeg{iSeg,1} = tst(begIdxCorr(iSeg):endIdxCorr(iSeg),:);
%     
% end


end

function [accProc, tstProc] = fillACC(acc, tst, ctrl, fs)
% pre-process acc data so that missing values are interpolated and
% inter-file spaces are NaNs spaced out according to common
% sampling-frequency of all files in session
% ctrl goes from 0 to N for a single file of N samples (some may be
% missing). For an NWB file where multiple files were concatenated, the 0
% to N ctrl count will reset and start again from 0 at the startint point
% of the next file. 


Ts = 1/fs;                           % Sampling Time
% tst = [0:10  14:20  25:30]*0.00274;       % Create Time Series With Gaps
% acc = rand(size(tst));                      % Create Matching Data
tst = tst';
tn = round(tst/Ts);                       % Create Indices With Gaps
dt = diff([0 tn]);                      % Create Vector Of Differences
tg = find(dt > 1);                      % Find Indices Of Gaps
gaps = dt(tg)-1;                        % Find Lengths Of Gaps
ti = linspace(min(tst),max(tst),max(tn));   % Create Continuous Time Vector
fi = interp1(tst,acc,ti);                   % Create Matching Data Vector
for k1 = 1:length(tg)                   % Loop Filling Gaps In �f� With NaN
    q = [tg(k1):tg(k1)+gaps(k1)-1];
    fi(tg(k1):tg(k1)+gaps(k1)-2) = nan(1,gaps(k1)-1);
end


% First subdivide 




nCtrls = length(ctrl);
for ictrl = 1:nCtrls
    
    
end

% then interpolate any short missing values within recordings


end

function [Sacc, f, tSpec] = spectrogramacc3(acc, fs, cfg)

% Detrend data, linear best fit is subtracted from each channel
acc = detrend(acc);


% Highpass filter out slow components

[b,a] = butter(3, cfg.fcDetrend/(fs/2), 'high');
accFilt = filtfilt(b, a, acc);



%% Get combined triaxial spectrogram

window = floor(cfg.TWIN * fs);

% Sum PSDs from all 3 accel channels
clear Z S
for ch = 1:3
    [Z(:,:,ch), f, tSpec] = spectrogram(accFilt(:,ch), window, cfg.NOVERLAP, [], fs, 'yaxis');
    % tRef = linspace(tSpec(1), tSpec(end), numel(tSpec));
    S(:,:,ch) = abs(Z(:,:,ch));
    S(:,:,ch) = S(:,:,ch).^2;

end
% Combine spectral power from all three axes (summate)
Sacc = sum(S, 3);

end

function ScleanCell = removeMovArtInAll(SaccCell, cfg)
% First concatenate all segments spectrograms together and keep track of
% their place within segment cells
nSegs = length(SaccCell);
SaccAll = [];
origIdx = [];
for iSeg = 1:nSegs
    iSacc = SaccCell{iSeg,1};
    origIdx = [origIdx, 1:size(iSacc, 2)]; % keeps tabls of time windows in spectrograms
    SaccAll = [SaccAll, SaccCell{iSeg,1}];   
    
end


% Identify movement artifact windows in spectrogram

[~, isMovArt] = util.remoutliers(sum(SaccAll, 1), ...
    'bound', 'upper', ...
    'MADthresh', cfg.MADthresh);
temp = SaccAll;
temp(:,isMovArt) = NaN;
SaccMarkOutliers = SaccAll;
SaccMarkOutliers(:,isMovArt) = NaN;


% Show Spectrogram with detected movement artifacts marked as 100's values
if cfg.displayWorkflowFigs
    f4 = figure;
    ax4 = axes;
    ax4.Parent = f4;
    surf(tSpec, f, 10*log10(SaccMarkOutliers), ...
        'Parent', ax4, ...
        'EdgeColor', 'none');
    hold(ax4, 'on');
    xlabel('Time (min)');
    grid on
    axis(ax4, 'tight')
    view(0,90)
    colorbar
    % ax4.YLim = [0, 300];
    % ax4.CLim = [-80, -30];
    title('Mov Artifacts marked');
    ylabel('Frequency (Hz)');
    f4.Position = [2555 57 560 421];

end

% Show Spectrogram with movement artifacts removed totally
if cfg.displayWorkflowFigs
    f4 = figure;
    ax4 = axes;
    ax4.Parent = f4;
%     tset = t;
%     tset(isMovArt) = [];
    surf(1:size(temp, 2), f, 10*log10(temp), ...
        'Parent', ax4, ...
        'EdgeColor', 'none');
    hold(ax4, 'on');
    xlabel('samples');
    grid on
    axis(ax4, 'tight')
    view(0,90)
    colorbar
    ax4.YLim = [0, 150];
    % ax4.CLim = [-80, -30];
    title('Mov Artifacts removed')
    ylabel('Frequency (Hz)');
    f4.Position = [2554 566 560 420];

end


% Get average PSD free of movement artifact windows
ScleanAll = temp; % movement artifact windows removed
% PSD = mean(ScleanAll, 2);
if cfg.displayWorkflowFigs
    f1 = figure; ax = axes;
    plot(f, (PSD));
    grid on; 
    ax.XLim = [0, 35];
    %     ax.YLim = [0, popMaxPSD];
    xlabel('Frequency (Hz)');
%     title(['Minutes since harmaline inject: ' num2str(minutes(sessTabRow.harRefTime(1)))], 'Interpreter', 'none')

end


% Output nan-marked spectrograms back in their original cells
begIdx = find(origIdx == 1);
endIdx = [(begIdx(2:end)-1), size(ScleanAll, 2)];
for iSeg = 1:nSegs
    ScleanCell{iSeg,1} = ScleanAll(:,(begIdx(iSeg):endIdx(iSeg)));
    
end

% Remove nan-marked columns in each spectrogram (pertaining to time-window
% labeled as movemend artifact) for each segment
for iSeg = 1:nSegs
    iSaccnan = ScleanCell{iSeg,1};
    nanIdx = isnan(sum(iSaccnan, 1));
    iSaccClean = iSaccnan;
    iSaccClean(:,nanIdx) = [];
    
    ScleanCell{iSeg,1} = iSaccClean;
    
end



end