% test script to derive lfp data from an epoch around voltage trace
clear
%%

SESSION_ID = 'TremoLfpDBS-191115-100127';
titStr = '8 mg/kg, 2nd';



procSuffix = '_LFP_C0-C2_';
BASELINE_SUFFIX = 'naivePre_exp';
WASHIN_SUFFIX = 'washin_exp';
% naiveSuffix = 'naivePre_exp';
nwbReadPn = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Processing\NWBdata\';

% number of samples for moving average in the "smooth" function to smooth
% the PSDs...
smoothSamps = 100;



%% Load data

% Get data from naive portion of current experiment day 
nwbBaseline = nwbRead([nwbReadPn SESSION_ID procSuffix BASELINE_SUFFIX '.nwb']);
lfpObjBaseline = nwbBaseline.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('lfp');

fsLfp = lfpObjBaseline.starting_time_rate;
dtLfp = 1/fsLfp;

dataLfpBaseline = nwbnrtl.util.loadTimeSeriesData(lfpObjBaseline);
tstLfpBaseline = nwbnrtl.util.loadTimeSeriesTimestamps(lfpObjBaseline);


% Get data from harmaline-washing part of current exp day
nwbWashin = nwbRead([nwbReadPn SESSION_ID procSuffix WASHIN_SUFFIX '.nwb']);
lfpObjWashin = nwbWashin.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('lfp');

dataLfpWashin = nwbnrtl.util.loadTimeSeriesData(lfpObjWashin);
tstLfpWashin = nwbnrtl.util.loadTimeSeriesTimestamps(lfpObjWashin);



%% Get PSD of naive LFPs averaged over 10-second windows

winTime = 10; % seconds
winSamp = floor(10 * fsLfp);

% Count how many non-overlapping 10-sec windows there are in naive LFPs
lenWashin = length(dataLfpBaseline);
winSampBinEdges = winSamp:winSamp:lenWashin;
winSampBins = [0, winSampBinEdges(1:(end-1))]';
winSampBins = winSampBins + 1;
winSampBins = [winSampBins, winSampBinEdges'];


dataWindow = dataLfpBaseline(winSampBins(1,1):winSampBins(1,2));

params.Fs = fsLfp;
params.tapers = [3, 5];


% get a quick sample of windowed data so we'll know the length
[S, f] = mtspectrumc(dataWindow, params);
nSampPSD = length(S); 

% Fill up a pre-allocated matrix to hold all the PSDs
nWindows = size(winSampBins, 1);
PSDs_base = zeros(nSampPSD, nWindows);
for iWin = 1:nWindows
    dataWindow = dataLfpBaseline(winSampBins(iWin,1):winSampBins(iWin,2));
    [PSDs_base(:,iWin), f] = mtspectrumc(dataWindow, params);
    
end

% get the average, and smooth it, to keep as baseline for later
% normalization
avPSDbaseline = mean(PSDs_base, 2);



% Remove spectral content around line noise and harmonics
remWin = [-0.3, 0.3]; % Hz window around central frequency to remove

lineNoiseHz = 60:60:(fsLfp/2);
isNaNIdx = zeros(1, length(S));
for iHz = lineNoiseHz
    % mark frequencies that fall within the removal window
    isRemiHz = (f > (remWin(1)+iHz)) & (f < (remWin(2)+iHz));
    isNaNIdx = isNaNIdx | isRemiHz;
    
end
isRemHigh = f > 400;
isRemLow = f < 0.5;
% isNaNIdx = isNaNIdx | isRemHigh | isRemLow;
% isNaNIdx = isNaNIdx | isRemHigh | isRemLow;
isRemFr = isRemHigh | isRemLow;
% fClean = f;
% SClean = S;
% fClean(isRemIdx) = NaN;

avPSDbaseline(isNaNIdx) = NaN;
avPSDbaseline(isRemFr) = [];
f(isRemFr) = [];

% SCleanDB = pow2db(SClean);
% SCleanDBsm = smooth(SCleanDB, 1000);
% SCleansm = smooth(SClean, 1000);
% figure; plot(fClean, SCleanDBsm)
% figure; plot(fClean, SCleanDB)

% figure; plot(f,avPSDbaseline)
% figure; plot(f, 10*log10(avPSDbaseline))
% figure; plot(f, 10*log10(smooth(avPSDbaseline, smoothSamps)));



%% Get PSDs of washin LFPs progressing in non-overlapping 10-second windows

winTime = 10; % seconds
winSamp = floor(10 * fsLfp);

% Count how many non-overlapping 10-sec windows there are in WASHIN LFPs
lenWashin = length(dataLfpWashin);
winSampBinEdges = winSamp:winSamp:lenWashin;
winSampBins = [0, winSampBinEdges(1:(end-1))]';
winSampBins = winSampBins + 1;
winSampBins = [winSampBins, winSampBinEdges'];


dataWindow = dataLfpWashin(winSampBins(1,1):winSampBins(1,2));

params.Fs = fsLfp;
params.tapers = [3, 5];


% get a quick sample of windowed data so we'll know the length
[S, f] = mtspectrumc(dataWindow, params);
nSampPSD = length(S); 

% Fill up a pre-allocated matrix to hold all the PSDs
nWindows = size(winSampBins, 1);
PSDs = zeros(nSampPSD, nWindows);
for iWin = 1:nWindows
    dataWindow = dataLfpWashin(winSampBins(iWin,1):winSampBins(iWin,2));
    [iPSD, f] = mtspectrumc(dataWindow, params);
    iPSD = smooth(iPSD, smoothSamps);
    PSDs(:,iWin) = iPSD;
    
end




% Remove spectral content around line noise and harmonics
remWin = [-0.3, 0.3]; % Hz window around central frequency to remove

lineNoiseHz = 60:60:(fsLfp/2);
isNaNIdx = zeros(1, length(S));
for iHz = lineNoiseHz
    % mark frequencies that fall within the removal window
    isRemiHz = (f > (remWin(1)+iHz)) & (f < (remWin(2)+iHz));
    isNaNIdx = isNaNIdx | isRemiHz;
    
end
isRemHigh = f > 400;
isRemLow = f < 0.5;
% isNaNIdx = isNaNIdx | isRemHigh | isRemLow;
% isNaNIdx = isNaNIdx | isRemHigh | isRemLow;
isRemFr = isRemHigh | isRemLow;
% fClean = f;
% SClean = S;
% fClean(isRemIdx) = NaN;

PSDs(isNaNIdx,:) = NaN;
PSDs(isRemFr,:) = [];
f(isRemFr) = [];


% % Smooth the PSDs and the average baseline
% PSDs = smoothdata(PSDs, 1, 'movmean', smoothSamps);
% avPSDbaseline = smoothdata(avPSDbaseline, 1, 'movmean', smoothSamps);
% figure; plot(f, avPSDbaseline);
% figure; plot(f, PSDs(:,1));

% Now scale each of these PSDs by the baseline average PSD, and take
% 10*log10 to get the decibel increase (or decrease) relative to naive
deltaPSDdb = PSDs;
for iWin = 1:nWindows
    deltaPSDdb(:,iWin) = 10*log10(PSDs(:,iWin) ./ avPSDbaseline);
    
end

% smooth the result 
deltaPSDdb = smoothdata(deltaPSDdb, 1, 'movmean', smoothSamps);

% Now plot each PSD on the same figure
f1 = figure; 
ax1 = axes;
hold on
grid on
ax1.XLim = [0, 50];
for iWin = 1:nWindows
    plot(f, deltaPSDdb(:,iWin));          
    
end
ylabel('PSD 10*log10(Pwashin / Pnaive)');
xlabel('Frequency (Hz)');
title([titStr ', 10-sec windows over time']);




%%
% 
% % get PSD of all washin data:
% 
% params.Fs = fsLfp;
% % fc = 100;
% % params.fpass = [0, 300];
% params.tapers = [3, 5];
% % movingwin = [5, 1];
% 
% 
% % [S, f] = mtspectrumc(lfp, params);
% % disp('begin spectrogram')
% % tic
% [S, f] = mtspectrumc(dataLfpBaseline, params);
% 
% % Remove spectral content around line noise and harmonics
% remWin = [-0.2, 0.2]; % Hz window around central frequency to remove
% 
% lineNoiseHz = 60:60:(fsLfp/2);
% isNaNIdx = zeros(1, length(S));
% for iHz = lineNoiseHz
%     % mark frequencies that fall within the removal window
%     isRemiHz = (f > (remWin(1)+iHz)) & (f < (remWin(2)+iHz));
%     isNaNIdx = isNaNIdx | isRemiHz;
%     
% end
% isRemHigh = f > 400;
% isRemLow = f < 0.5;
% isNaNIdx = isNaNIdx | isRemHigh | isRemLow;
% fClean = f;
% SClean = S;
% fClean(isNaNIdx) = [];
% SClean(isNaNIdx) = [];
% 
% SCleanDB = pow2db(SClean);
% SCleanDBsm = smooth(SCleanDB, 1000);
% SCleansm = smooth(SClean, 1000);
% figure; plot(fClean, SCleanDBsm)
% figure; plot(fClean, SCleanDB)
% 
% % ax = gca; ax.XScale = 'log'
% 
% 
% %% Fit a model to the line to use for flattening the PSD
% 
% [fitResult, gof] = createCurveFit(fClean, SCleanDBsm)
% 
% % Smooth the PSD of interest
% 
% % Now Get the flattened PSD by dividing normal PSD by model
% 
% 
% 
% 
% f1 = figure; 
% plot(f, pow2db(S))
% 
% % run spectrogram on it
% window = 2^10;
% noverlap = floor(0.25*window);
% 
% % spectrogram for naivePre portion
% [Z, f, t] = spectrogram(dataLfpTot, window, noverlap, [], fsLfp, 'yaxis');
% tRef = linspace(tstLfpTot(1), tstLfpTot(end), numel(t));
% fset = f(f < 80);
% Zset = Z(1:length(fset),:);
% S = abs(Zset);
% ax = axes;
% ax.Parent = f1;
% surf(tRef/60, fset, pow2db(S.^2), ...
%     'Parent', ax, ...
%     'EdgeColor', 'none');
% hold(ax, 'on');
% xlabel('Time (min)');
% grid on
% axis(ax, 'tight')
% view(0,90)
% 
% 
% % 
% % % spectrogram for naivePre portion
% % [Z, f, t] = spectrogram(dataLfpNaive, window, noverlap, [], fsLfp, 'yaxis');
% % S = abs(Z);
% % ax(1) = subplot(1,2,1);
% % ax(1).Parent = f1;
% % surf(t/60, f, pow2db(S.^2), ...
% %     'Parent', ax(1), ...
% %     'EdgeColor', 'none');
% % hold(ax(1), 'on');
% % xlabel('Time (min)');
% % grid on
% % axis(ax(1), 'tight')
% % view(0,90)
% 
% % 
% % % spectrogram for washin portion
% % [Z, f, t] = spectrogram(dataLfpWashin, window, noverlap, [], fsLfp, 'yaxis');
% % S = abs(Z);
% % ax(2) = subplot(1,2,2);
% % ax(2).Parent = f1;
% % surf(t/60, f, pow2db(S.^2), ...
% %     'Parent', ax(2), ...
% %     'EdgeColor', 'none');
% % hold(ax(2), 'on');
% % xlabel('Time (min)');
% % grid on
% % axis(ax(1), 'tight')
% % view(0,90)
% % 
% % linkaxes(ax)
% set(ax, 'YLim', [0, 80]);
% set(ax, 'YTick', 0:10:150);
% set(ax, 'CLimMode', 'manual');
% set(ax, 'CLim', [-75 -60]);
% 
% % grid(ax(:),'on');
% % axis(ax(:),'tight');
% % ax(:).YLim = [0, 80];
% % ax(:).YTick = 0:10:150;
% 
% % 
% % % ax(:).CLimMode = 'manual';
% % ax(:).CLim = [-180 -20]; % [-156.5356 -77.2316]
% % ax(:).CLim = [-75 -60]; % [-156.5356 -77.2316]
% 
% colorbar
% title(titStr)
% 
% 
% %% 
% 
% % Look at PSD's over time
% t = (1/fsLfp) * (0:(length(dataLfpNaive)-1));
% % secWindow = floor(300 * fsLfp)
% % plot psd from 1 to 300 seconds
% % [pxx, f] = pwelch(dataLfp(1:secWindow), fsLfp);
% f2 = figure;
% ax2 = axes;
% % pwelch(dataLfp(1:secWindow), fsLfp);
% pwin = 2^14;
% poverlap = pwin/2;
% pnfft = pwin;
% [pxx, f] = pwelch(dataLfpNaive, pwin, poverlap, pnfft, fsLfp);
% plot(f, 10*log10(pxx))
% ylabel('PSD (dB/Hz)')
% xlabel('Frequency (Hz)')
% ax2.XScale = 'log'
% ax2.XLim = [1, 150];
% grid minor
% 
% 
% %% generate raw nwb file
