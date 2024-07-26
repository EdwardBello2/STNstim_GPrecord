% test script to derive lfp data from an epoch around voltage trace
clear; close all
%%

SESSION_ID = 'TremoLfpDBS-191011-104322';
dose = 6;
expos = 1;
titStr = [num2str(dose) ' mg/kg, expos ' num2str(expos)];



procSuffix = '_LFP_C0-C2_';
BASELINE_SUFFIX = 'naivePre_exp';
% WASHIN_SUFFIX = 'washin_exp';
% naiveSuffix = 'naivePre_exp';
nwbReadPn = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Processing\NWBdata\';



% Load data

% Get data from naive portion of current experiment day 
nwbBaseline = nwbRead([nwbReadPn SESSION_ID procSuffix BASELINE_SUFFIX '.nwb']);
lfpObjBaseline = nwbBaseline.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('lfp');

fsLfp = lfpObjBaseline.starting_time_rate;
dtLfp = 1/fsLfp;

dataLfpBaseline = nwbnrtl.util.loadTimeSeriesData(lfpObjBaseline);
tstLfpBaseline = nwbnrtl.util.loadTimeSeriesTimestamps(lfpObjBaseline);


% Specify parameters for spectrogram:
tWin = 5; % seconds
window = floor(tWin * fsLfp);
overlapPerc = 25;
noverlap = floor((overlapPerc/100) * window);


% Prep save information:
savePn = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Processing\naiveExp_LFPs\';
saveFn = [procSuffix(2:end-1), '_dose', num2str(dose), ...
    '_expos', num2str(expos), ...
     '_', BASELINE_SUFFIX, '_specgram', ...
     '_win', num2str(tWin), 's', ...
     '_overlap', num2str(overlapPerc), 'perc'];

 
 
% Get spectrogram of entire of naive portion, and average PSD from it

% Good ol' fashioned FFT matlab default method of spectrogram
f1 = figure; 

% spectrogram for naivePre portion
[Z, f, t] = spectrogram(dataLfpBaseline, window, noverlap, [], fsLfp, 'yaxis');
tRef = linspace(t(1), t(end), numel(t));
% fset = f(f < 80);
% Zset = Z(1:length(fset),:);
Zset = Z;
Sbase = abs(Zset);
Sbase = Sbase.^2;
ax = axes;
ax.Parent = f1;
surf(t/60, f, 10*log10(Sbase), ...
    'Parent', ax, ...
    'EdgeColor', 'none');
hold(ax, 'on');
xlabel('Time (min)');
grid on
axis(ax, 'tight')
view(0,90)
colorbar
ax.YLim = [0, 300];
ax.CLim = [-80, -30];
title([titStr ', naivePre Spectrogram'])
ylabel('Frequency (Hz)');
f1.Position = [1990 567 560 420];
% Get average power (S) PSD over time

%% Fit

% get subset of frequency content from spectrogram
outlierFreqRange = [1, 300]; % Hz
isFset = (f >= outlierFreqRange(1)) & (f < outlierFreqRange(2));
SbaseSet = (Sbase(isFset,:)); % rescaling in order to put different frequency amplitudes closer to each other


[~, isOutlier] = util.remoutliers(sum(SbaseSet, 1), 'bound', 'upper');
temp = Sbase;
temp(:,isOutlier) = [];
SbaseMarkOutliers = Sbase;
SbaseMarkOutliers(:,isOutlier) = 100;

f4 = figure;
ax4 = axes;
ax4.Parent = f4;
surf(t/60, f, 10*log10(SbaseMarkOutliers), ...
    'Parent', ax4, ...
    'EdgeColor', 'none');
hold(ax4, 'on');
xlabel('Time (min)');
grid on
axis(ax4, 'tight')
view(0,90)
colorbar
ax4.YLim = [0, 300];
ax4.CLim = [-80, -30];
title(['outliers marked'])
ylabel('Frequency (Hz)');
f4.Position = [2554 566 560 420];


% 
% % figure; histogram(zscore(sum(Sbase, 1)))
% % temp = Sbase;
% 
% % first remove obvious mov artifact outliers
% isRem = zscore(sum(Sbase, 1)) > 2; % mark any PSDs with power greater than 2 stdvs 
% isOutlier = isOutlier | isRem;
% temp(:,isRem) = [];
% figure; histogram(zscore(sum(temp, 1)));
% 
% % next remove outliers strictly by exceeding 3 stdv
% isRem = zscore(sum(temp, 1)) > 3; % mark any PSDs with power greater than 2 stdvs 
% isOutlier = isOutlier | isRem;
% figure; histogram(zscore(sum(temp, 1)));
CIperc = 95; % percent
PSDbaseAv = mean(temp, 2);
% get accompanying Confidence Interval view of PSD 
nFreqs = length(PSDbaseAv);
CI = zeros(nFreqs, 2);
scaleTemp = 10*log10(temp);
for iFr = 1:nFreqs
    CI(iFr,1:2) = util.getConfInt(scaleTemp(iFr,:), CIperc);
    
end


f2 = figure; ax = axes; 
plot(f, 10*log10(PSDbaseAv));
hold on; plot(f, (CI(:,1)), 'r'); plot(f, (CI(:,2)), 'r')
ax.XScale = 'log'; grid on
title([titStr ', clean mean PSD (' num2str(CIperc) '% CI)'])
f2.Position = [1989 59 1133 420];
ax.YLim = [-100, -20];



%% Normalize naive spectrogram by its own average PSD, display as surface

% Normalize
PSDbaseNorm = Sbase ./ (repmat(PSDbaseAv, 1, size(Sbase, 2)));

sm = 500;
f3 = figure; 
ax = axes;
ax.Parent = f3;
surf(tRef/60, f, 10*log10(PSDbaseNorm), ...
    'Parent', ax, ...
    'EdgeColor', 'none');
hold(ax, 'on');
xlabel('Time (min)');
grid on
axis(ax, 'tight')
view(0,90)

colorbar
ax.YLim = [0, 300];
ax.CLim = [-9, 9];
title([titStr ', naivePre Gain-spectrogram'])
f3.Position = [3116 566 560 420];



%% Save results for later

% save the:
% 1) naive spectrogram (with time and freq)
% 2) TF array for detected outliers PSDs


save([savePn saveFn], 'Sbase', 't', 'f', 'isOutlier')





