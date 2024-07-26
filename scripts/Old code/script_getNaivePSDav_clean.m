% test script to derive lfp data from an epoch around voltage trace
clear; close all
%% Display Spectrograms of the naive data

SESSION_ID = 'TremoLfpDBS-191115-100127';
dose = 8;
expos = 2;
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

%% Get Truly Clean naive PSD for normalizing later washin spectrogram (gain-o-gram)

% % Remove movement artifacts
% outlierFreqRange = [1, 1000]; % Hz
% isFset = (f >= outlierFreqRange(1)) & (f < outlierFreqRange(2));
% SbaseSet = (Sbase(isFset,:)); % rescaling in order to put different frequency amplitudes closer to each other

[~, isMovArt] = util.remoutliers(sum(Sbase, 1), ...
    'bound', 'upper', ...
    'MADthresh', 2);
temp = Sbase;
temp(:,isMovArt) = [];
SbaseMarkOutliers = Sbase;
SbaseMarkOutliers(:,isMovArt) = 100;

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
title(['Mov Artifacts marked'])
ylabel('Frequency (Hz)');
f4.Position = [2555 57 560 421];

f4 = figure;
ax4 = axes;
ax4.Parent = f4;
tset = t;
tset(isMovArt) = [];
surf([1:size(temp, 2)], f, 10*log10(temp), ...
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
title(['Mov Artifacts removed'])
ylabel('Frequency (Hz)');
f4.Position = [2554 566 560 420];




% Fit line to Average PSD of naive epoch (free of mov artifacts)
fitFreqRange = [1, 300];
isFset = (f >= fitFreqRange(1)) & (f < fitFreqRange(2));
tempPSDav = mean(temp, 2);
tempPSDav(~isFset) = [];
fcFit = f(isFset);


% Get flattened naive PSD model
a = 0.00014;
b = 0;
c = 1;
d = 0;
n = 1.4;
y = modelPSD(fcFit, n, a, b, c, d); % custom in script
figure; ax = axes;
plot(fcFit, 10*log10(tempPSDav)); hold on; plot(fcFit, 10*log10(y), 'r');
ax.XScale = 'log';
title('naive PSD av, Movement artifacts removed')
fig = gcf; fig.Position = [1992 59 560 420];
flatPSD = tempPSDav ./ y;
figure; plot(fcFit, 10*log10(flatPSD))


% Get flattend PSD for spectrogram frequency regions of interest 
outlierFreqRange = [2, 45]; % Hz
isFset = (f >= outlierFreqRange(1)) & (f < outlierFreqRange(2));
tempFset = temp(isFset,:);
fset = f(isFset);
modelPSDFset = modelPSD(fset, n, a, b, c, d); % custom in script

naiveSpecgramNorm = (tempFset ./ (repmat(modelPSDFset, 1, size(tempFset, 2))));

% % Detect transient power outliers now (FINALLY)
% [~, isOutlier] = util.remoutliers(sum(naiveSpecgramNorm, 1), ...
%     'bound', 'upper', ...
%     'MADthresh', 2);
% temp = naiveSpecgramNorm;
% temp(:,isOutlier) = [];
% figure; histogram(sum(naiveSpecgramNorm, 1))
% % SbaseMarkOutliers = naiveSpecgramNorm;
% % SbaseMarkOutliers(:,isOutlier) = 100;




% Alternative way to tell when band power is outlier: do outlier detection
% individually for each and every freq band over time
nFreqs = size(naiveSpecgramNorm, 1);
isOutlier = false(nFreqs, size(naiveSpecgramNorm, 2));
for iFreq = 1:nFreqs
    [~,isOutlier(iFreq,:)] = util.remoutliers(naiveSpecgramNorm(iFreq,:), ...
    'bound', 'both', ...
    'MADthresh', 6);

end

isBadTime = any(isOutlier,1);


% Finish by getting truly clean average PSD for normalizing later washin
noArts = Sbase;
noArts(:,isMovArt) = [];
noOutliers = noArts;
noOutliers(:,isBadTime) = [];

cleanPSDav = mean(noOutliers, 2);
figure; plot(f, 10*log10(cleanPSDav));
ax = gca; ax.XScale = 'log';
title('Clean naive PSD average')


% % get subset of frequency content from spectrogram
% outlierFreqRange = [1, 1000]; % Hz
% isFset = (f >= outlierFreqRange(1)) & (f < outlierFreqRange(2));
% SbaseSet = (Sbase(isFset,:)); % rescaling in order to put different frequency amplitudes closer to each other
% 
% 
% [~, isOutlier] = util.remoutliers(sum(SbaseSet, 1), 'bound', 'upper');
% temp = Sbase;
% temp(:,isOutlier) = [];
% SbaseMarkOutliers = Sbase;
% SbaseMarkOutliers(:,isOutlier) = 100;
% 
% f4 = figure;
% ax4 = axes;
% ax4.Parent = f4;
% surf(t/60, f, 10*log10(SbaseMarkOutliers), ...
%     'Parent', ax4, ...
%     'EdgeColor', 'none');
% hold(ax4, 'on');
% xlabel('Time (min)');
% grid on
% axis(ax4, 'tight')
% view(0,90)
% colorbar
% ax4.YLim = [0, 300];
% ax4.CLim = [-80, -30];
% title(['outliers marked'])
% ylabel('Frequency (Hz)');
% f4.Position = [2554 566 560 420];


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
% CIperc = 95; % percent
% PSDbaseAv = mean(temp, 2);
% % get accompanying Confidence Interval view of PSD 
% nFreqs = length(PSDbaseAv);
% CI = zeros(nFreqs, 2);
% scaleTemp = 10*log10(temp);
% for iFr = 1:nFreqs
%     CI(iFr,1:2) = util.getConfInt(scaleTemp(iFr,:), CIperc);
%     
% end
% 
% 
% f2 = figure; ax = axes; 
% plot(f, 10*log10(PSDbaseAv));
% hold on; plot(f, (CI(:,1)), 'r'); plot(f, (CI(:,2)), 'r')
% ax.XScale = 'log'; grid on
% title([titStr ', clean mean PSD (' num2str(CIperc) '% CI)'])
% f2.Position = [1989 59 1133 420];
% ax.YLim = [-100, -20];



%% Normalize naive spectrogram by its own average PSD, display as surface

% Normalize
SbaseClean = Sbase;
SbaseClean(:,isMovArt) = [];
SbaseClean(:,isBadTime) = [];
nTimeWins = size(SbaseClean, 2)

PSDbaseNorm = SbaseClean ./ (repmat(cleanPSDav, 1, nTimeWins));

f3 = figure; 
ax = axes;
ax.Parent = f3;
surf([1:nTimeWins], f, 10*log10(PSDbaseNorm), ...
    'Parent', ax, ...
    'EdgeColor', 'none');
hold(ax, 'on');
xlabel('window num');
grid on
axis(ax, 'tight')
view(0,90)

colorbar
ax.YLim = [0, 80];
ax.CLim = [-9, 9];
title([titStr ', naivePre Gain-spectrogram'])
f3.Position = [2552 568 560 420];

% 
% f4 = figure;
% ax4 = axes;
% ax4.Parent = f4;
% tset = t;
% tset(isMovArt) = [];
% surf([1:size(temp, 2)], f, 10*log10(temp), ...
%     'Parent', ax4, ...
%     'EdgeColor', 'none');
% hold(ax4, 'on');
% xlabel('Time (min)');
% grid on
% axis(ax4, 'tight')
% view(0,90)
% colorbar
% ax4.YLim = [0, 300];
% ax4.CLim = [-80, -30];
% title(['Mov Artifacts removed'])
% ylabel('Frequency (Hz)');
% f4.Position = [2554 566 560 420];



%% Save results for later

% save the:
% 1) naive spectrogram (with time and freq)
% 2) TF array for detected outliers PSDs


save([savePn saveFn], 'Sbase', 't', 'f', 'isMovArt', 'isBadTime')




function y = modelPSD(x, n, a, b, c, d)


% Get flattened naive PSD

y = (a./((c*(x+b)).^n)) + d;

end


