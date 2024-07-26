% test script to derive lfp data from an epoch around voltage trace
clear; close all
%% Get the right data

SESSION_ID = 'TremoLfpDBS-191115-100127';
dose = 8;
expos = 2;
chPair = 'C4-C6'; % str, 'C0-C2' | 'C2-C4' | 'C4-C6'

% % 
% % % load naive-baseline PSDaverage for normalization purposes
projRootPn = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';
% naivePSDnormPn = 'Data Processing\naiveExp_LFPs\';
% T = readtable([projRootPn naivePSDnormPn 'naiveExp_LFPs_metadata.xlsx']);
% % 


procSuffix = ['_LFP_' chPair '_'];
BASELINE_SUFFIX = 'naivePre_exp';
WASHIN_SUFFIX = 'washin_exp';
% naiveSuffix = 'naivePre_exp';
nwbReadPn = [projRootPn 'Data Processing\NWBdata\'];

% number of samples for moving average in the "smooth" function to smooth
% the PSDs...
smoothSamps = 100;
titStr = [procSuffix ', ' num2str(dose) ' mg/kg, expos ' num2str(expos)];

reportSubDir = 'Reports\Report-201019_doseStudyLFPs_multi\';


% % Get naive PSD average
% 
% iRow = (T.dose_mgkg_ == dose) & (T.exposure == expos);
% fn = T.filename{iRow};
% load([projRootPn naivePSDnormPn fn], 'Sbase', 't', 'f', 'isMovArt', 'isMovArt', 'isBadTime');
% 
% % Remove outlier time windows from spectrogram (Sbase) and get mean PSD
% Sbase(:,isMovArt) = [];
% Sbase(:,isBadTime) = [];
% PSDbaseAv = mean(Sbase, 2);
% % figure; ax = axes;
% % plot(f, 10*log10(PSDbaseAv));
% % ax.XScale = 'log';
% 


% Load washin data

% iRow = (T.dose_mgkg_ == dose) & (T.exposure == expos);
% nwbfn = T.ancestorFileName{iRow};
% nwbpn = T.ancestorFilePath{iRow};
% Get data from harmaline-washing part of current exp day
nwbWashin = nwbRead([nwbReadPn SESSION_ID procSuffix WASHIN_SUFFIX '.nwb']);
lfpObjWashin = nwbWashin.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('lfp');

dataLfpWashin = nwbnrtl.util.loadTimeSeriesData(lfpObjWashin);
tstLfpWashin = nwbnrtl.util.loadTimeSeriesTimestamps(lfpObjWashin);

fsLfp = lfpObjWashin.starting_time_rate;
dtLfp = 1/fsLfp;

% run spectrogram on it
tWin = 5; % seconds
window = floor(tWin * fsLfp);
noverlap = floor(0.25*window);



% Get spectrogram of washin portion, normalized to naive portion, and display as surface


% Good ol' fashioned FFT matlab default method of spectrogram

% % run spectrogram on it
% tWin = 5; % seconds
% window = floor(tWin * fsLfp);
% noverlap = floor(0.25*window);

% spectrogram for naivePre portion
[Z, f, t] = spectrogram(dataLfpWashin, window, noverlap, [], fsLfp, 'yaxis');
tRef = linspace(t(1), t(end), numel(t));
% fset = f(f < 80);
% Zset = Z(1:length(fset),:);
Zset = Z;
Swash = abs(Zset);
Swash = Swash.^2;

% f1 = figure; 
% ax = axes;
% ax.Parent = f1;
% surf(tRef/60, f, pow2db(Swash), ...
%     'Parent', ax, ...
%     'EdgeColor', 'none');
% hold(ax, 'on');
% xlabel('Time (min)');
% grid on
% axis(ax, 'tight')
% view(0,90)
% colorbar
% ax.YLim = [0, 300];


% Get average power (S) PSD over time
% PSDWashNorm = Swash ./ (repmat(PSDbaseAv, 1, size(Swash, 2)));
PSDWashNorm = Swash;

sm = 500;
f2 = figure; 
ax = axes;
ax.Parent = f2;
surf(tRef/60, f, 10*log10(PSDWashNorm), ...
    'Parent', ax, ...
    'EdgeColor', 'none');
hold(ax, 'on');
xlabel('Time since har inj (min)');
grid on
axis(ax, 'tight')
view(0,90)
colorbar
ax.CLim = [-80, -35];
ax.YLim = [0, 100];
a1 = gca;
title([titStr, ', Washin, CLim: [' num2str(a1.CLim(1)) ' ' num2str(a1.CLim(2)) ']'], 'interpreter', 'none')
ylabel('Frequency (Hz)')
set(gcf, 'Position', [2096 189 1390 726])

% Save this figure in the report folder
savefig(gcf, [projRootPn, reportSubDir, 'WashinSpecgram' ...
    '_', chPair, ...
    '_dose', num2str(dose), ...
    '_expos', num2str(expos)]);

