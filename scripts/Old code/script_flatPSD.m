% test script to derive lfp data from an epoch around voltage trace
clear
%%

SESSION_ID = 'TremoLfpDBS-191115-100127';
titStr = '8 mg/kg, 2nd';

% tdtDataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';
% TDTtank = 'Uva-191115';
% TDTblockFirst = 'TremoLfpDBS-191115-100127';
procSuffix = '_LFP_C0-C2_';
epochSuffix = 'naivePre_exp';
% naiveSuffix = 'naivePre_exp';
nwbReadPn = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Processing\NWBdata\';



% read NWB file for current experiment day
% nwbNiave = nwbRead([nwbReadPn SESSION_ID nwbSuffix naiveSuffix '.nwb']);
% lfpObjNaive = nwbNiave.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('lfp');

nwbEpoch = nwbRead([nwbReadPn SESSION_ID procSuffix epochSuffix '.nwb']);
lfpObjEpoch = nwbEpoch.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('lfp');


fsLfp = lfpObjEpoch.starting_time_rate;
dtLfp = 1/fsLfp;

% dataLfpNaive = nwbnrtl.util.loadTimeSeriesData(lfpObjNaive);
% tstLfpNaive = nwbnrtl.util.loadTimeSeriesTimestamps(lfpObjNaive);

dataLfpEpoch = nwbnrtl.util.loadTimeSeriesData(lfpObjEpoch);
tstLfpEpoch = nwbnrtl.util.loadTimeSeriesTimestamps(lfpObjEpoch);

% % re-reference timestamps and combine them
% tstLfpWashin = tstLfpWashin - tstLfpWashin(1);
% tstLfpNaive = tstLfpNaive - tstLfpNaive(end) - dtLfp;
% 
% dataLfpTot = [dataLfpNaive; dataLfpWashin];
% tstLfpTot = [tstLfpNaive; tstLfpWashin];
% 

% get PSD of all washin data:

params.Fs = fsLfp;
% fc = 100;
% params.fpass = [0, 300];
params.tapers = [3, 5];
% movingwin = [5, 1];


% [S, f] = mtspectrumc(lfp, params);
% disp('begin spectrogram')
% tic
[S, f] = mtspectrumc(dataLfpEpoch, params);

% Remove spectral content around line noise and harmonics
remWin = [-0.2, 0.2]; % Hz window around central frequency to remove

lineNoiseHz = 60:60:(fsLfp/2);
isRemIdx = zeros(1, length(S));
for iHz = lineNoiseHz
    % mark frequencies that fall within the removal window
    isRemiHz = (f > (remWin(1)+iHz)) & (f < (remWin(2)+iHz));
    isRemIdx = isRemIdx | isRemiHz;
    
end
isRemHigh = f > 400;
isRemLow = f < 0.5;
isRemIdx = isRemIdx | isRemHigh | isRemLow;
fClean = f;
SClean = S;
fClean(isRemIdx) = [];
SClean(isRemIdx) = [];

SCleanDB = pow2db(SClean);
SCleanDBsm = smooth(SCleanDB, 1000);
SCleansm = smooth(SClean, 1000);
figure; plot(fClean, SCleanDBsm)
figure; plot(fClean, SCleanDB)

% ax = gca; ax.XScale = 'log'


%% Fit a model to the line to use for flattening the PSD

[fitResult, gof] = createCurveFit(fClean, SCleanDBsm)

% Smooth the PSD of interest

% Now Get the flattened PSD by dividing normal PSD by model




f1 = figure; 
plot(f, pow2db(S))

% run spectrogram on it
window = 2^10;
noverlap = floor(0.25*window);

% spectrogram for naivePre portion
[Z, f, t] = spectrogram(dataLfpTot, window, noverlap, [], fsLfp, 'yaxis');
tRef = linspace(tstLfpTot(1), tstLfpTot(end), numel(t));
fset = f(f < 80);
Zset = Z(1:length(fset),:);
S = abs(Zset);
ax = axes;
ax.Parent = f1;
surf(tRef/60, fset, pow2db(S.^2), ...
    'Parent', ax, ...
    'EdgeColor', 'none');
hold(ax, 'on');
xlabel('Time (min)');
grid on
axis(ax, 'tight')
view(0,90)


% 
% % spectrogram for naivePre portion
% [Z, f, t] = spectrogram(dataLfpNaive, window, noverlap, [], fsLfp, 'yaxis');
% S = abs(Z);
% ax(1) = subplot(1,2,1);
% ax(1).Parent = f1;
% surf(t/60, f, pow2db(S.^2), ...
%     'Parent', ax(1), ...
%     'EdgeColor', 'none');
% hold(ax(1), 'on');
% xlabel('Time (min)');
% grid on
% axis(ax(1), 'tight')
% view(0,90)

% 
% % spectrogram for washin portion
% [Z, f, t] = spectrogram(dataLfpWashin, window, noverlap, [], fsLfp, 'yaxis');
% S = abs(Z);
% ax(2) = subplot(1,2,2);
% ax(2).Parent = f1;
% surf(t/60, f, pow2db(S.^2), ...
%     'Parent', ax(2), ...
%     'EdgeColor', 'none');
% hold(ax(2), 'on');
% xlabel('Time (min)');
% grid on
% axis(ax(1), 'tight')
% view(0,90)
% 
% linkaxes(ax)
set(ax, 'YLim', [0, 80]);
set(ax, 'YTick', 0:10:150);
set(ax, 'CLimMode', 'manual');
set(ax, 'CLim', [-75 -60]);

% grid(ax(:),'on');
% axis(ax(:),'tight');
% ax(:).YLim = [0, 80];
% ax(:).YTick = 0:10:150;

% 
% % ax(:).CLimMode = 'manual';
% ax(:).CLim = [-180 -20]; % [-156.5356 -77.2316]
% ax(:).CLim = [-75 -60]; % [-156.5356 -77.2316]

colorbar
title(titStr)


%% 

% Look at PSD's over time
t = (1/fsLfp) * (0:(length(dataLfpNaive)-1));
% secWindow = floor(300 * fsLfp)
% plot psd from 1 to 300 seconds
% [pxx, f] = pwelch(dataLfp(1:secWindow), fsLfp);
f2 = figure;
ax2 = axes;
% pwelch(dataLfp(1:secWindow), fsLfp);
pwin = 2^14;
poverlap = pwin/2;
pnfft = pwin;
[pxx, f] = pwelch(dataLfpNaive, pwin, poverlap, pnfft, fsLfp);
plot(f, 10*log10(pxx))
ylabel('PSD (dB/Hz)')
xlabel('Frequency (Hz)')
ax2.XScale = 'log'
ax2.XLim = [1, 150];
grid minor


%% generate raw nwb file
