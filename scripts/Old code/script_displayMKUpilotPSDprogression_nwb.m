% script for displaying ACC PSDs over the data-segments in time for each
% subject. Purpose: to visually see the PSD shape change. As for the
% time-element, that's covered in script_displayMKUpilotTremorInTime_nwb.m,
% focusing on the tremor-power/non-tremor power ratio over time. 

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

errbCapSize = 0;
errbLineWidth = 1;
errbC = 'k';

scatSz = 30;
scatC = 'k'; 
scatCsig = 'r';

harLineWidth = 1.5;

pchipLineStyle = '--';
linC_M = 'b';
linC_K = 'm';
linC_U = 'g';



%% Mela's data

% Mela Pilot data of choice: 12 mg/kg day
nwbFullPath = [nwbDataPath 'MelaPilot\'];
% SESSION_ID = 'harPilotMela-150928-094900'; % 2 mg/kg day
% SESSION_ID = 'harPilotMela-150930-082000'; % 6 mg/kg day
SESSION_ID = 'harPilotMela-151013-072100'; % 12 mg/kg day

% Load nwb master file for this session, get acc TimeSeries object
disp('Reading NWB file...')
nwbM = nwbRead([nwbFullPath SESSION_ID '_MASTER.nwb']);
disp('DONE!')

% Calculate spectrogram-based PSD windows from segments of recorded data
[psdSegsM, tstSegHarRefM, fSegsM, resultsInfo] = func_gatherACCproc_allPSDTimeSegs_har(nwbM, cfg);

nSegs = length(psdSegsM);
for iSeg = 1:nSegs
    psdDispM{iSeg,1} = mean(psdSegsM{iSeg}, 2, 'omitnan');
    tmidM(iSeg,1) = (tstSegHarRefM{iSeg,1}(1) + tstSegHarRefM{iSeg,1}(end)) / 2; 

end
tmidM = tmidM / 60; 

% Get frequency peaks for each data segment
fspan = [0 35];
[psdFpksSegsM, psdPpksSegsM] = getSegmentFpeaks(psdSegsM, fSegsM, fspan);


% Display frequency peak summary statistic with confidence intervals for
% each segment in time
[~,FpksAvM,cinegM,ciposM,~,~] = func_prepACCproc_fbandTimeSegs_errorbar(...
    psdFpksSegsM, tstSegHarRefM, cfg, 'mean');


s = 1;
figure('Position', [1985 61 560 420]); 
axFpks(s) = axes;
% ax(s) = subplot(3,1,s);
errorbar(tmidM,FpksAvM,cinegM,ciposM,[],[], '.',  ...
    'CapSize', errbCapSize, 'LineWidth', errbLineWidth,'Color', errbC);
grid on; hold on;
xqM = tmidM(1):0.1:tmidM(end);
yyM = pchip(tmidM, FpksAvM, xqM);
pM = plot(xqM, yyM, 'LineStyle', pchipLineStyle, 'Color', linC_M);
scatter(tmidM, FpksAvM, scatSz, scatC, 'filled');
ylabel('Frequency (Hz)')
title('Mela frequency peaks')



% Display peak PSD power summary statistic with confidence intervals for
% each segment in time
[~,PpksAvM,cinegM,ciposM,~,~] = func_prepACCproc_fbandTimeSegs_errorbar(...
    psdPpksSegsM, tstSegHarRefM, cfg, 'median');

% Normalizing step for power measures
normVal = max(PpksAvM);
PpksAvM = PpksAvM / normVal;
cinegM = cinegM / normVal;
ciposM = ciposM / normVal;

s = 1;
figure('Position', [1984 565 560 420]); 
axPpks(s) = axes;
% ax(s) = subplot(3,1,s);
errorbar(tmidM,PpksAvM,cinegM,ciposM,[],[], '.',  ...
    'CapSize', errbCapSize, 'LineWidth', errbLineWidth,'Color', errbC);
grid on; hold on;
xqM = tmidM(1):0.1:tmidM(end);
yyM = pchip(tmidM, PpksAvM, xqM);
pM = plot(xqM, yyM, 'LineStyle', pchipLineStyle, 'Color', linC_M);
scatter(tmidM, PpksAvM, scatSz, scatC, 'filled');
ylabel('Normalized peak power')
title('Mela power peaks')



%% Kiwi's data

% Kiwi Pilot data of choice: 10 mg/kg day
nwbFullPath = [nwbDataPath 'KiwiPilot\'];
SESSION_ID = 'harPilotKiwi-160517-122121';

% Load nwb master file for this session, get acc TimeSeries object
disp('Reading NWB file...')
nwbK = nwbRead([nwbFullPath SESSION_ID '_MASTER.nwb']);
disp('DONE!')

% Calculate spectrogram-based PSD windows from segments of recorded data
[psdSegsK, tstSegHarRefK, fSegsK, resultsInfo] = func_gatherACCproc_allPSDTimeSegs_har(nwbK, cfg);

nSegs = length(psdSegsK);
for iSeg = 1:nSegs
    psdDispK{iSeg,1} = mean(psdSegsK{iSeg}, 2, 'omitnan');
    tmidK(iSeg,1) = (tstSegHarRefK{iSeg,1}(1) + tstSegHarRefK{iSeg,1}(end)) / 2; 
    
end
tmidK = tmidK/60; % minutes


% Get frequency peaks for each data segment
fspan = [10 35];
[psdFpksSegsKhi, psdPpksSegsKhi] = getSegmentFpeaks(psdSegsK, fSegsK, fspan);


% Display frequency peak summary statistic with confidence intervals for
% each segment in time
[~,FpksAvKhi,cinegKhi,ciposKhi,tnegKhi,tposKhi] = func_prepACCproc_fbandTimeSegs_errorbar(...
    psdFpksSegsKhi, tstSegHarRefK, cfg, 'mean');


% Get frequency peaks for each data segment
fspan = [0 10];
[psdFpksSegsKlo, psdPpksSegsKlo] = getSegmentFpeaks(psdSegsK, fSegsK, fspan);


% Display frequency peak summary statistic with confidence intervals for
% each segment in time
[~,FpksAvKlo,cinegKlo,ciposKlo,tnegKlo,tposKlo] = func_prepACCproc_fbandTimeSegs_errorbar(...
    psdFpksSegsKlo, tstSegHarRefK, cfg, 'mean');


s = 2;
figure('Position', [2549 62 560 420]); 
axFpks(s) = axes;
% ax(s) = subplot(3,1,s);
grid on; hold on;


% Hi freq tremor
errorbar(tmidK,FpksAvKhi,cinegKhi,ciposKhi,[],[], '.',  ...
    'CapSize', errbCapSize, 'LineWidth', errbLineWidth,'Color', errbC);
xqKhi = tmidK(1):0.1:tmidK(end);
yyKhi = pchip(tmidK, FpksAvKhi, xqKhi);
pKhi = plot(xqKhi, yyKhi, 'LineStyle', pchipLineStyle, 'Color', linC_K);
scatter(tmidK, FpksAvKhi, scatSz, scatC, 'filled');


% Lo freq tremor
errorbar(tmidK,FpksAvKlo,cinegKlo,ciposKlo,[],[], '.',  ...
    'CapSize', errbCapSize, 'LineWidth', errbLineWidth,'Color', errbC);
xqKlo = tmidK(1):0.1:tmidK(end);
yyKlo = pchip(tmidK, FpksAvKlo, xqKlo);
pKlo = plot(xqKlo, yyKlo, 'LineStyle', pchipLineStyle, 'Color', linC_K);
scatter(tmidK, FpksAvKlo, scatSz, scatC, 'filled');


ylabel('Frequency (Hz)')
title('Kiwi frequency peaks')


figure('Position', [2547 567 560 420]); 

% Display peak PSD power summary statistic with confidence intervals for
% each segment in time
% Hi
[~,PpksAvKhi,cinegKhi,ciposKhi,tnegKhi,tposKhi] = func_prepACCproc_fbandTimeSegs_errorbar(...
    psdPpksSegsKhi, tstSegHarRefK, cfg, 'median');

% Normalizing step for power measures
normVal = max(PpksAvKhi);
PpksAvKhi = PpksAvKhi / normVal;
cinegKhi = cinegKhi / normVal;
ciposKhi = ciposKhi / normVal;

s = 2;
axPpks(s) = subplot(2,1,1);
% ax(s) = subplot(3,1,s);
errorbar(tmidK,PpksAvKhi,cinegKhi,ciposKhi,[],[], '.',  ...
    'CapSize', errbCapSize, 'LineWidth', errbLineWidth,'Color', errbC);
grid on; hold on;
xqKhi = tmidK(1):0.1:tmidK(end);
yyKhi = pchip(tmidK, PpksAvKhi, xqKhi);
pKhi = plot(xqKhi, yyKhi, 'LineStyle', pchipLineStyle, 'Color', linC_K);
scatter(tmidK, PpksAvKhi, scatSz, scatC, 'filled');

% Lo
[~,PpksAvKlo,cinegKlo,ciposKlo,tnegKlo,tposKlo] = func_prepACCproc_fbandTimeSegs_errorbar(...
    psdPpksSegsKlo, tstSegHarRefK, cfg, 'median');

% Normalizing step for power measures
normVal = max(PpksAvKlo);
PpksAvKlo = PpksAvKlo / normVal;
cinegKlo = cinegKlo / normVal;
ciposKlo = ciposKlo / normVal;

s = 4;
axPpks(s) = subplot(2,1,2);
% ax(s) = subplot(3,1,s);
errorbar(tmidK,PpksAvKlo,cinegKlo,ciposKlo,[],[], '.',  ...
    'CapSize', errbCapSize, 'LineWidth', errbLineWidth,'Color', errbC);
grid on; hold on;
xqKlo = tmidK(1):0.1:tmidK(end);
yyKlo = pchip(tmidK, PpksAvKlo, xqKlo);
pKlo = plot(xqKlo, yyKlo, 'LineStyle', pchipLineStyle, 'Color', linC_K);
scatter(tmidK, PpksAvKlo, scatSz, scatC, 'filled');



ylabel('Normalized peak power')
title('Kiwi power peaks')



%% Uva's data

% Uva Pilot data of choice: 8 mg/kg day
nwbFullPath = [nwbDataPath 'UvaPilot\'];
% SESSION_ID = 'harPilotUva-180219-132618'; % 8 mg/kg day
SESSION_ID = 'harPilotUva-180221-113120'; % 12 mg/kg day
% SESSION_ID = 'harPilotUva-180305-103016'; % 12 mg/kg day

% Load nwb master file for this session, get acc TimeSeries object
disp('Reading NWB file...')
nwbU = nwbRead([nwbFullPath SESSION_ID '_MASTER.nwb']);
disp('DONE!')

% Calculate spectrogram-based PSD windows from segments of recorded data
[psdSegsU, tstSegHarRefU, fSegsU, resultsInfo] = func_gatherACCproc_allPSDTimeSegs_har(nwbU, cfg);

nSegs = length(psdSegsU);
for iSeg = 1:nSegs
    psdDispU{iSeg,1} = mean(psdSegsU{iSeg}, 2, 'omitnan');
    tmidU(iSeg,1) = (tstSegHarRefU{iSeg,1}(1) + tstSegHarRefU{iSeg,1}(end)) / 2; 
    
end
tmidU = tmidU / 60;

% Get frequency peaks for each data segment
fspan = [0 35];
[psdFpksSegsU, psdPpksSegsU] = getSegmentFpeaks(psdSegsU, fSegsU, fspan);


% Display frequency peak summary statistic with confidence intervals for
% each segment in time
[~,FpksAvU,cinegU,ciposU,tnegU,tposU] = func_prepACCproc_fbandTimeSegs_errorbar(...
    psdFpksSegsU, tstSegHarRefU, cfg, 'mean');


s = 3;
figure('Position', [3112 63 560 420]); 
axFpks(s) = axes;
% ax(s) = subplot(3,1,s);
errorbar(tmidU,FpksAvU,cinegU,ciposU,[],[], '.',  ...
    'CapSize', errbCapSize, 'LineWidth', errbLineWidth,'Color', errbC);
grid on; hold on;

xqU = tmidU(1):0.1:tmidU(end);
yyU = pchip(tmidU, FpksAvU, xqU);
pU = plot(xqU, yyU, 'LineStyle', pchipLineStyle, 'Color', linC_U);
scatter(tmidU, FpksAvU, scatSz, scatC, 'filled');
ylabel('Frequency (Hz)')
title('Uva frequency peaks')


% Display peak PSD power summary statistic with confidence intervals for
% each segment in time
[~,PpksAvU,cinegU,ciposU,tnegU,tposU] = func_prepACCproc_fbandTimeSegs_errorbar(...
    psdPpksSegsU, tstSegHarRefU, cfg, 'median');

% Normalizing step for power measures
normVal = max(PpksAvU);
PpksAvU = PpksAvU / normVal;
cinegU = cinegU / normVal;
ciposU = ciposU / normVal;

s = 3;
figure('Position', [3110 568 560 420]); 
axPpks(s) = axes;
% ax(s) = subplot(3,1,s);
errorbar(tmidU,PpksAvU,cinegU,ciposU,[],[], '.',  ...
    'CapSize', errbCapSize, 'LineWidth', errbLineWidth,'Color', errbC);
grid on; hold on;

xqU = tmidU(1):0.1:tmidU(end);
yyU = pchip(tmidU, PpksAvU, xqU);
pU = plot(xqU, yyU, 'LineStyle', pchipLineStyle, 'Color', linC_U);
scatter(tmidU, PpksAvU, scatSz, scatC, 'filled');
ylabel('Normalized peak power')
title('Uva power peaks')



%% Display All results
% RGB gray: [128 128 128]
% RGB dark-orange: [255 140 0]

% view3d = [18.6000 36.2000];
view3d = [16.2000 73];
freqLim = [0 35];
linWid = 1.0;
sm = 1;


fig1 = figure('Position', [1986 569 1717 420]); 

% Mela
s = 1;
% fig(s) = figure('Position', [1990 420 560 420]); 
axPSD(s) = subplot(1,3,s);
hold on
nSegs = length(psdDispM);
% cM = flipud(parula(nSegs));
cM = parula(nSegs);
for iSeg = 1:nSegs
    psd = psdDispM{iSeg};
%     psd = smooth(psdDispM{iSeg}, 5);
    plot3(fSegsM{iSeg}, ones(length(psd), 1) * (tmidM(iSeg)), smooth(psdDispM{iSeg}, sm), ...
        'LineWidth', linWid, 'Color', cM(iSeg,:));
    
end
xlabel('Frequency (Hz)')
% ylabel('Time (minutes)')
title('S1 (Mela)')


% Kiwi
s = 2;
% fig(s) = figure('Position', [2553 420 560 420]); 
axPSD(s) = subplot(1,3,s);
hold on
nSegs = length(psdDispK);
% cK = flipud(parula(nSegs));
cK = parula(nSegs);
for iSeg = 1:nSegs
    psd = psdDispK{iSeg};
%     psd = smooth(psdDispK{iSeg}, 5);
    plot3(fSegsK{iSeg}, ones(length(psd), 1) * (tmidK(iSeg)), smooth(psdDispK{iSeg}, sm), ...
        'LineWidth', linWid, 'Color', cK(iSeg,:));
    
end
xlabel('Frequency (Hz)')
% ylabel('Time (minutes)')
title('S2 (Kiwi)')


% Uva
s = 3;
% fig(s) = figure('Position', [3115 420 560 420]); 
axPSD(s) = subplot(1,3,s);
hold on
nSegs = length(psdDispU);
% cU = flipud(parula(nSegs));
cU = parula(nSegs);
for iSeg = 1:nSegs
    psd = psdDispU{iSeg};
%     psd = smooth(psdDispU{iSeg}, 5);
    plot3(fSegsU{iSeg}, ones(length(psd), 1) * (tmidU(iSeg)), smooth(psdDispU{iSeg}, sm), ...
        'LineWidth', linWid, 'Color', cU(iSeg,:));
    
end
xlabel('Frequency (Hz)')
% ylabel('Time (minutes)')
title('S3 (Uva)')


set(axPSD(:), 'XLim', freqLim, 'View', view3d, ...
    'ZTickLabel', [], 'ZTick', [], ...
    'YLim', [-10 150])


set(axFpks(:), 'XLim', [-10 150], 'YLim', [0 22])
for i = 1:length(axFpks)
    plot(axFpks(i), [0 0], [0 22], 'r', 'LineWidth', harLineWidth);
    
end


% for i = 1:length(axPpks)
%     plot(axPpks(i), [0 0], [0 22], 'r', 'LineWidth', harLineWidth);
%     
% end
set(axPpks(:), 'XLim', [-10 150], 'YLim', [0 1.5]);
% set(axPpks(:), 'XLim', [-10 150]);



%% SUB-FUNCTIONS

function [psdFpksSegsM, psdPpksSegsM] = getSegmentFpeaks(psdSegsM, fSegsM, fspan)

% Get frequency peaks for each data segment
% fspan = [0 35];
nSegs = length(psdSegsM);
for iSeg = 1:nSegs
    % Specify subselection of frequencies to look at
    SegF = fSegsM{iSeg,1};
    isRange = (SegF >= fspan(1)) & (SegF < fspan(2));
    SegFsub = SegF(isRange);
    
    
    % Get frequency peaks
    SegPSDs = psdSegsM{iSeg,1}(isRange,:);
    nPSDs = size(SegPSDs, 2);
    for iPSD = 1:nPSDs
        Ppk(iPSD) = max(SegPSDs(:,iPSD));
        Fpk(iPSD) = SegFsub(SegPSDs(:,iPSD) == Ppk(iPSD));
        
    end
    
    psdPpksSegsM{iSeg,1} = Ppk;
    psdFpksSegsM{iSeg,1} = Fpk;
    
    
end



end

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