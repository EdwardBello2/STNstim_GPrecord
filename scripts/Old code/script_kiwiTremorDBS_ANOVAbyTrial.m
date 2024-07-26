% script to statistically test if tremors changed during DBS stim in Kiwi's data, based on
% "script_displayKiwiStimTremorEffect.m"

% Collected on microstrain G-link hardware
clear; close all

% INPUTS

projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';

% for loading
dataAcqPn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\';
directoryPn = 'Kiwi\KiwiExperiment_05_17_2016\';
inputMetaData = 'RecordingsMetadata_selectDisplay.mat';

% for saving
SavePn = 'Data Processing\Kiwi_tremors\';


% CONSTANTS

% Downsampling high sampling frequency files from ~8kHz down to ~512 Hz, to
% be near in sampling frequency to the 500 Hz sampling frequency files:
DWNFACTOR = 16;

% Spectrogram parameters:
TWIN = 0.5; % seconds, sliding time-window of spectrogram
NOVERLAP = 0; % samples, no overlaps for these spectrogam windows



%% Load metadata table for all recordings 

% pre-load the metadata table with file info; use only the "DBS" triaType
load([projRootPath 'Data Acquisition\Kiwi Experiment\' inputMetaData], 'Tselect');

isDBS = strcmp('DBS',Tselect.trialType);
Tselect(~isDBS,:) = [];



%% Loop thru all individual files adding values to table

nRows = height(Tselect);
tic
for iRow = 1:nRows
    close all
    % Load in data and get into format of "acc" and "fs"
    fn = Tselect.filename{iRow};
    load([projRootPath 'Data Acquisition\Kiwi Experiment\' fn], ...
        'data', 'time', 'sampling_frequency', 'output_signal');


    % if it's a stim recording, the first channel will be a recording of TACS
    % or DBS signal, for a total of 4 channels; otherwise there should just be
    % 3 channels of accel
    if size(data, 2) == 4 
        acc = data(:,2:4);

    else
        acc = data(:,1:3);

    end
    t = time;
    fs = sampling_frequency;

    % Detrend data, linear best fit is subtracted from each channel
    acc = detrend(acc);


    % Highpass filter out slow components
    fc = 1; % Hz
    [b,a] = butter(3, fc/(fs/2), 'high');
    accFilt = filtfilt(b, a, acc);

    % Because of stupidity I recorded things at two different sampling
    % frequencies. Now, PSD peak amplitudes are significantly affected by 
    % sampling frequency, so to properly compare data recorded with
    % different fs we need to downsample the faster one to get as close to
    % the slower as possible. 

    % two frequencies used: 500 Hz, 8,193.4 Hz
    if fs > 500
    %     DWNFACTOR = 16;
        fs = fs/DWNFACTOR; % bring down to 512.0852
        t = downsample(t, DWNFACTOR);
        clear y
        for i = 1:3
            y(:,i) = decimate(accFilt(:,i), DWNFACTOR);

        end
        accFilt = y;

    else
    %     dTime = t;
    end

    % figure; plot(t, accFilt)
    % title('highpassed')



    %% Get combined triaxial spectrogram

    window = floor(TWIN * fs);

    % Sum PSDs from all 3 accel channels
    clear Z S
    for ch = 1:3
        [Z(:,:,ch), f, tSpec] = spectrogram(accFilt(:,ch), window, NOVERLAP, [], fs, 'yaxis');
        % tRef = linspace(tSpec(1), tSpec(end), numel(tSpec));
        S(:,:,ch) = abs(Z(:,:,ch));
        S(:,:,ch) = S(:,:,ch).^2;

    end
    % Combine spectral power from all three axes (summate)
    Sacc = sum(S, 3);

%     % Get average power (S) PSD over time
%     % PSDWashNorm = Swash ./ (repmat(PSDbaseAv, 1, size(Swash, 2)));
%     f2 = figure; 
%     ax = axes;
%     ax.Parent = f2;
%     surf(tSpec, f, 10*log10(Sacc), ...
%         'Parent', ax, ...
%         'EdgeColor', 'none');
%     hold(ax, 'on');
%     % xlabel('Time (min)');
%     grid on
%     axis(ax, 'tight')
%     view(0,90)
%     colorbar
%     % ax.CLim = [-9, 9];
%     % ax.YLim = [0, 20];
%     % title([titStr, ', Washin, dBgain (Washin/Naive)'])



    %% Identify movement artifact windows in spectrogram

    [~, isMovArt] = util.remoutliers(sum(Sacc, 1), ...
        'bound', 'upper', ...
        'MADthresh', 3);
    temp = Sacc;
    temp(:,isMovArt) = [];
    SaccMarkOutliers = Sacc;
    SaccMarkOutliers(:,isMovArt) = 100;
    tSpecClean = tSpec; tSpecClean(isMovArt) = [];
    Sclean = temp;


%     % Show Spectrogram with detected movement artifacts marked as 100's values
%     f4 = figure;
%     ax4 = axes;
%     ax4.Parent = f4;
%     surf(tSpec, f, 10*log10(SaccMarkOutliers), ...
%         'Parent', ax4, ...
%         'EdgeColor', 'none');
%     hold(ax4, 'on');
%     xlabel('Time (min)');
%     grid on
%     axis(ax4, 'tight')
%     view(0,90)
%     colorbar
%     % ax4.YLim = [0, 300];
%     % ax4.CLim = [-80, -30];
%     title('Mov Artifacts marked');
%     ylabel('Frequency (Hz)');
%     f4.Position = [2555 57 560 421];
%     
%     
%     % Show Spectrogram with movement artifacts removed totally
%     f4 = figure;
%     ax4 = axes;
%     ax4.Parent = f4;
%     tset = t;
%     tset(isMovArt) = [];
%     surf(1:size(temp, 2), f, 10*log10(temp), ...
%         'Parent', ax4, ...
%         'EdgeColor', 'none');
%     hold(ax4, 'on');
%     xlabel('samples');
%     grid on
%     axis(ax4, 'tight')
%     view(0,90)
%     colorbar
%     ax4.YLim = [0, 150];
%     % ax4.CLim = [-80, -30];
%     title('Mov Artifacts removed')
%     ylabel('Frequency (Hz)');
%     f4.Position = [2554 566 560 420];
% 
% 
%     % Get average PSD free of movement artifact windows
%     Sclean = temp; % movement artifact windows removed
%     PSD = mean(Sclean, 2);
%     
%     f1 = figure; ax = axes;
%     plot(f, (PSD));
%     grid on; 
%     ax.XLim = [0, 35];
% %     ax.YLim = [0, popMaxPSD];
%     xlabel('Frequency (Hz)');
%     title(['Minutes since harmaline inject: ' num2str(minutes(Tselect.harRefTime(iRow)))], 'Interpreter', 'none')



%     %% Split time-series data into two bands
% 
%     % binsPerDecade = 2;
%     % nBins = 40;
%     % binEdgeLeft = 10^(0);
%     % binEdges = binEdgeLeft*10.^((0:nBins)/binsPerDecade)
% 
%     fcLo = [3.1623, 10]; % Hz
%     fcHi = [10, 31.6228]; % Hz 
% 
%     % Low band: 
%     [b, a] = butter(2, fcLo / (fs/2), 'bandpass');
%     accLo = filtfilt(b, a, accFilt);
% 
%     % High band:
%     [b, a] = butter(2, fcHi / (fs/2), 'bandpass');
%     accHi = filtfilt(b, a, accFilt);
% 
% 
%     %% Get PSD estimates for both bands as above
% 
%     % LOW-BAND
%     window = floor(TWIN * fs);
%     clear Z S
%     for ch = 1:3
%         [Z(:,:,ch), f, tSpec] = spectrogram(accLo(:,ch), window, NOVERLAP, [], fs, 'yaxis');
%         % tRef = linspace(tSpec(1), tSpec(end), numel(tSpec));
%         S(:,:,ch) = abs(Z(:,:,ch));
%         S(:,:,ch) = S(:,:,ch).^2;
% 
%     end
%     SaccLo = sum(S, 3);
%     SaccLoClean = SaccLo;
%     SaccLoClean(:,isMovArt) = [];
%     PSDlo = mean(SaccLoClean, 2);
%     % figure; plot(f, 10*log10(PSDlo));
% 
% 
%     % HIGH-BAND
%     window = floor(TWIN * fs);
%     clear Z S
%     for ch = 1:3
%         [Z(:,:,ch), f, tSpec] = spectrogram(accHi(:,ch), window, NOVERLAP, [], fs, 'yaxis');
%         % tRef = linspace(tSpec(1), tSpec(end), numel(tSpec));
%         S(:,:,ch) = abs(Z(:,:,ch));
%         S(:,:,ch) = S(:,:,ch).^2;
% 
%     end
%     SaccHi = sum(S, 3);
%     SaccHiClean = SaccHi;
%     SaccHiClean(:,isMovArt) = [];
%     PSDhi = mean(SaccHiClean, 2);
%     % figure; plot(f, 10*log10(PSDhi));
%     tSpec = tSpec';
% 
%     % TRACK peak frequency results for table
%     fPeakLo(iRow,1) = f(PSDlo == max(PSDlo));
%     fPeakHi(iRow,1) = f(PSDhi == max(PSDhi));
%     
%     % TRACK PSD peak amplitudes for each band too
%     aPeakLo(iRow,1) = max(PSDlo);
%     aPeakHi(iRow,1) = max(PSDhi);

% 
% 
%     %% Get hilbert instantaneous band power for lo and hi
% 
%     accLoMag = abs(hilbert(accLo));
%     % accLoPha = angle(hilbert(accLo));
%     accHiMag = abs(hilbert(accHi)); 
%     % accHiPha = angle(hilbert(accHi)); 
% 
% 
% 
%     %% Remove data points pertaining to movment artifact time windows, detected
%     % from the Spectrogram work...
% 
%     accLoMagSum = sum(accLoMag, 2);
%     accHiMagSum = sum(accHiMag, 2);
% 
%     % Go from known mov artifact specgram time windows to labeling time-series
%     % indices to exclude from final analysis:
% 
%     % specify m x n matrix to track which time-series samples to mark. m:
%     % samples in time-series; n: number of spectrogram windows detected as mov
%     % artifact
% 
%     nArts = sum(isMovArt);
%     idxArt = find(isMovArt);
%     artWinIdx = false(length(accLoMag), sum(isMovArt));
% 
%     for iArt = 1:nArts
%         % specgram window time represents center of window
%         i_timeArt = tSpec(idxArt(iArt));
%         i_timeRange = [(i_timeArt - TWIN/2), (i_timeArt + TWIN/2)];
% 
%         % the values of time-series time that fall within this range will be
%         % marked as artifact sample indices
%         artWinIdx(:,iArt) = (t >= i_timeRange(1)) & (t < i_timeRange(2));
% 
%     end
%     % Collapse all detected mov artifact indices to 1D to mark the samples in
%     % time-series that need to be ignored
%     isMovArtTimeSeries = any(artWinIdx, 2);
% 
%     % figure; plot(accLoMagSum); hold on; plot(isMovArtTimeSeries);
%     % figure; plot(accHiMagSum); hold on; plot(isMovArtTimeSeries);
% 
%     % Remove values that fall within movement artifact time
%     accLoMagSum(isMovArtTimeSeries) = [];
%     accHiMagSum(isMovArtTimeSeries) = [];
%     t(isMovArtTimeSeries) = [];
% 
% 
% 
%     %% Get average Power for each band, and TRACK
% 
%     % Mean Power for total time
%     avLoMag(iRow,1) = mean(accLoMagSum);
%     avHiMag(iRow,1) = mean(accHiMagSum);
% 
%     % Cross-correlation for 
%     x = accLoMagSum;
%     y = accHiMagSum;
% 
%     [r, p] = corrcoef(x, y);
%     R(iRow,1) = r(2);
% 
%     x_boot = x;
% 
% 
%     % get p-test based on this empirical bootstrapp'd distribution
%     nBoots = 1000;
%     nSamps = length(x);
%     r_boot = zeros(nBoots, 1);
%     for iBoot = 1:nBoots
%         x_boot = circshift(x, randi(nSamps));
%         rtemp = corrcoef(x_boot, y);
%         r_boot(iBoot) = rtemp(2);
% 
%     end
%     Pc(iRow,1) = sum(r_boot > R(iRow,1)) / numel(r_boot);
% 
%     % figure; plot(x); hold on; plot(y); 
%     % title(['Pearson R: ' num2str(R) ', p = ' num2str(Pc)]); 
%     % legend('low band', 'high band')
    
    
    
    %% plot out instantaneous magnitude, and show portions of time where DBS was on
    
    tremorBand = [3, 20];
    isBand = (f >= tremorBand(1)) & (f < tremorBand(2));
    
    % specify parts of spectrogram time that belong to each epoch
    isPre = (tSpecClean >= 0) & (tSpecClean < 15);
    isDBS = (tSpecClean >= 15) & (tSpecClean < 45);
    isPos = (tSpecClean >= 45);
    
    % Get total bandpower for each time-window of spectrogram 
    prePSDs = Sclean(:,isPre);
    dbsPSDs = Sclean(:,isDBS);
    posPSDs = Sclean(:,isPos);
 
    % store each collection of points in appropriate cell (table row) 
    preTremor{iRow,1} = sum(prePSDs(isBand,:), 1);
    dbsTremor{iRow,1} = sum(dbsPSDs(isBand,:), 1);
    posTremor{iRow,1} = sum(posPSDs(isBand,:), 1);
    
    
%     
%     accLoMagSumBase = mean(accLoMagSum(isPre));
%     std1 = std(accLoMagSum(isPre));
%     accHiMagSumBase = mean(accHiMagSum(isPre));
%     std2 = std(accHiMagSum(isPre));
%     
%     
%     % plot instantaneous band power over time (non-normalized)
%     f1 = figure; ax = axes;
%     plot(t, accLoMagSum, '--r', t, accHiMagSum, '--b', 'LineWidth', 0.25); legend('Low band', 'High band');
%     ymax = ax.YLim(2);
%     hold on; plot([15, 15], [0, ymax], 'k', 'LineWidth', 2);
%     plot([45, 45], [0, ymax], 'k', 'LineWidth', 2);
%     ax.YLim = [0, ymax];
%     xlabel('time (seconds)'); ylabel('A.U.');
%     % show smoothed version too
%     span = 2500;
%     plot(t, smooth(accLoMagSum, span), 'r', t, smooth(accHiMagSum, span), 'b', 'LineWidth', 1)
%     accPreRawMean(iRow,1) = mean(accLoMagSum(isPre));
%     accPreRawStdv(iRow,1) = std(accLoMagSum(isPre));
%     accDbsRawMean(iRow,1) = mean(accLoMagSum(isDBS));
%     accDbsRawStdv(iRow,1) = std(accLoMagSum(isDBS));
%     accPosRawMean(iRow,1) = mean(accLoMagSum(isPos));
%     accPosRawStdv(iRow,1) = std(accLoMagSum(isPos));
%     title(fn, 'interpreter', 'none')
% %     title([fn, ' (av/std), pre:(', num2str(accPre(1), 2), '/', num2str(accPre(2), 2), ...
% %         '), dbs:(', num2str(accDBS(1), 2), '/', num2str(accDBS(2), 2), ...
% %         '), post:(', num2str(accPos(1), 2), '/', num2str(accPos(2), 2), ')'], ...
% %         'interpreter', 'none')
%     f1.Position = [117 495 560 420];
%     
% 
%     % plot band power over time normalized to pre-DBS baseline (%)
%     f2 = figure; ax = axes;
%     accLoMagSumNorm = 100*(accLoMagSum/avLoMagSessionBaseline);
%     accHiMagSumNorm = 100*(accHiMagSum/avHiMagSessionBaseline);
%     plot(t, accLoMagSumNorm, '--r', ...
%         t, accHiMagSumNorm, '--b', 'LineWidth', 0.25); 
%     grid on
%     ymin = ax.YLim(1); % stupid lims keep changing with each plot you add...
%     ymax = ax.YLim(2);
%     hold on; 
%     plot([15, 15], [ymin, ymax], 'k', 'LineWidth', 2);
%     plot([45, 45], [ymin, ymax], 'k', 'LineWidth', 2);
%     ax.YLim = [ymin, ymax];
%     xlabel('time (seconds)'); ylabel('Change from baseline (%)');
%     legend('Low band', 'High band');
%     % show smoothed version too
%     span = 2500;
%     plot(t, smooth(accLoMagSumNorm, span), 'r', ...
%     t, smooth(accHiMagSumNorm, span), 'b', 'LineWidth', 1);
%     accPrePercMean(iRow,1) = mean(accLoMagSumNorm(isPre));
%     accPrePercStdv(iRow,1) = std(accLoMagSumNorm(isPre));
%     accDbsPercMean(iRow,1) = mean(accLoMagSumNorm(isDBS));
%     accDbsPercStdv(iRow,1) = std(accLoMagSumNorm(isDBS));
%     accPosPercMean(iRow,1) = mean(accLoMagSumNorm(isPos));
%     accPosPercStdv(iRow,1) = std(accLoMagSumNorm(isPos));
%     title(fn, 'interpreter', 'none')
% %     title([fn, ' (av/std), pre:(', num2str(accPre(1), 2), '/', num2str(accPre(2), 2), ...
% %         '), dbs:(', num2str(accDBS(1), 2), '/', num2str(accDBS(2), 2), ...
% %         '), post:(', num2str(accPos(1), 2), '/', num2str(accPos(2), 2), ')'], ...
% %         'interpreter', 'none')
%     f2.Position = [679 494 560 420];
% 
% 
%     % plot band power over time normalized to pre-DBS baseline (dB)
%     f3 = figure; ax = axes;
%     accLoMagSumNorm = 10*log10(accLoMagSum/avLoMagSessionBaseline);
%     accHiMagSumNorm = 10*log10(accHiMagSum/avHiMagSessionBaseline);
%     plot(t, accLoMagSumNorm, '--r', ...
%         t, accHiMagSumNorm, '--b', 'LineWidth', 0.25); 
%     grid on
%     ymin = ax.YLim(1); % stupid lims keep changing with each plot you add...
%     ymax = ax.YLim(2);
%     hold on; 
%     plot([15, 15], [ymin, ymax], 'k', 'LineWidth', 2);
%     plot([45, 45], [ymin, ymax], 'k', 'LineWidth', 2);
%     ax.YLim = [ymin, ymax];
%     xlabel('time (seconds)'); ylabel('Change from baseline (dB)');
%     legend('Low band', 'High band');
%     % show smoothed version too
%     span = 2500;
%     plot(t, smooth(accLoMagSumNorm, span), 'r', ...
%     t, smooth(accHiMagSumNorm, span), 'b', 'LineWidth', 1);
%     accPreDbMean(iRow,1) = mean(accLoMagSumNorm(isPre));
%     accPreDbStdv(iRow,1) = std(accLoMagSumNorm(isPre));
%     accDbsDbMean(iRow,1) = mean(accLoMagSumNorm(isDBS));
%     accDBSDbStdv(iRow,1) = std(accLoMagSumNorm(isDBS));
%     accPosDbMean(iRow,1) = mean(accLoMagSumNorm(isPos));
%     accPosDbStdv(iRow,1) = std(accLoMagSumNorm(isPos));
%     title(fn, 'interpreter', 'none')
% %     title([fn, ' (av/std), pre:(', num2str(accPre(1), 2), '/', num2str(accPre(2), 2), ...
% %         '), dbs:(', num2str(accDBS(1), 2), '/', num2str(accDBS(2), 2), ...
% %         '), post:(', num2str(accPos(1), 2), '/', num2str(accPos(2), 2), ')'], ...
% %         'interpreter', 'none')
%     f3.Position = [1241 493 560 420];
    
    

end
toc



%% Get baseline average tremor power to normalize all tremor data 

% Load in data and get into format of "acc" and "fs"
fn = 'tremor_frequency_measure_05-17-16_time1220_arm.mat';
% fn = 'tremor_frequency_measure_05-17-16_time1233_arm.mat';
load([projRootPath 'Data Acquisition\Kiwi Experiment\' fn], ...
    'data', 'time', 'sampling_frequency', 'output_signal');


% if it's a stim recording, the first channel will be a recording of TACS
% or DBS signal, for a total of 4 channels; otherwise there should just be
% 3 channels of accel
if size(data, 2) == 4 
    acc = data(:,2:4);

else
    acc = data(:,1:3);

end
t = time;
fs = sampling_frequency;

% Detrend data, linear best fit is subtracted from each channel
acc = detrend(acc);


% Highpass filter out slow components
fc = 1; % Hz
[b,a] = butter(3, fc/(fs/2), 'high');
accFilt = filtfilt(b, a, acc);

% Because of stupidity I recorded things at two different sampling
% frequencies. Now, PSD peak amplitudes are significantly affected by 
% sampling frequency, so to properly compare data recorded with
% different fs we need to downsample the faster one to get as close to
% the slower as possible. 

% two frequencies used: 500 Hz, 8,193.4 Hz
if fs > 500
%     DWNFACTOR = 16;
    fs = fs/DWNFACTOR; % bring down to 512.0852
    t = downsample(t, DWNFACTOR);
    clear y
    for i = 1:3
        y(:,i) = decimate(accFilt(:,i), DWNFACTOR);

    end
    accFilt = y;

else
%     dTime = t;
end


% Get combined triaxial spectrogram

window = floor(TWIN * fs);

% Sum PSDs from all 3 accel channels
clear Z S
for ch = 1:3
    [Z(:,:,ch), f, tSpec] = spectrogram(accFilt(:,ch), window, NOVERLAP, [], fs, 'yaxis');
    % tRef = linspace(tSpec(1), tSpec(end), numel(tSpec));
    S(:,:,ch) = abs(Z(:,:,ch));
    S(:,:,ch) = S(:,:,ch).^2;

end
% Combine spectral power from all three axes (summate)
Sacc = sum(S, 3);


% Identify movement artifact windows in spectrogram

[~, isMovArt] = util.remoutliers(sum(Sacc, 1), ...
    'bound', 'upper', ...
    'MADthresh', 3);
temp = Sacc;
temp(:,isMovArt) = [];
SaccMarkOutliers = Sacc;
SaccMarkOutliers(:,isMovArt) = 100;
tSpecClean = tSpec; tSpecClean(isMovArt) = [];


% plot out instantaneous magnitude, and show portions of time where DBS was on

tremorBand = [3, 20];
isBand = (f >= tremorBand(1)) & (f < tremorBand(2));

baselineTremor = mean(sum(Sclean(isBand,:)))

% % specify parts of spectrogram time that belong to each epoch
% isPre = (tSpecClean >= 0) & (tSpecClean < 15);
% isDBS = (tSpecClean >= 15) & (tSpecClean < 45);
% isPos = (tSpecClean >= 45);
% 
% % Get total bandpower for each time-window of spectrogram 
% prePSDs = Sclean(:,isPre);
% dbsPSDs = Sclean(:,isDBS);
% posPSDs = Sclean(:,isPos);
% 
% % store each collection of points in appropriate cell (table row) 
% preTremor{iRow,1} = sum(prePSDs(isBand,:), 1);
% dbsTremor{iRow,1} = sum(dbsPSDs(isBand,:), 1);
% posTremor{iRow,1} = sum(posPSDs(isBand,:), 1);



%%

% add in collected data to table
Tdisp = Tselect;
Tdisp = [Tdisp, table(preTremor), table(dbsTremor), table(posTremor)];



% Choose table subselection for each statistical test, normalize
for iRow = 1:height(Tdisp)
    close all
    dataPre = 100 * Tdisp.preTremor{iRow} / baselineTremor;
    dataDbs = 100 * Tdisp.dbsTremor{iRow} / baselineTremor;
    dataPos = 100 * Tdisp.posTremor{iRow} / baselineTremor;
    
    % Perform ANOVA to test difference
    preLab = cell(numel(dataPre), 1);
    preLab(:) = {'pre'};
    dbsLab = cell(numel(dataDbs), 1);
    dbsLab(:) = {'dbs'};
    posLab = cell(numel(dataPos), 1);
    posLab(:) = {'pos'};
    labs = [preLab; dbsLab; posLab];
    
    dataAll = [dataPre'; dataDbs'; dataPos'];
    
    [p, tbl, stats] = anova1(dataAll, labs)
    
    
    % display barplots of means for each condition
    means(1) = mean(dataPre);
    means(2) = mean(dataDbs);
    means(3) = mean(dataPos);
    c = categorical({'Pre-DBS', 'DBS-on', 'Post-DBS'});
    c = reordercats(c, {'Pre-DBS', 'DBS-on', 'Post-DBS'});
    figure; br = bar(c, means);
    xlabel('trial epoch');
    ylabel('Tremor power (%, norm to baseline)');
    title(Tdisp.filename{iRow}, 'interpreter', 'none')
    
    figHandles = findobj('Type', 'figure');
    figHandles(1).Position = [3061         527         560         420];
    figHandles(2).Position = [2498         500         560         420];
    figHandles(3).Position = [1935         500         560         420];
    
    
end


% % perform ANOVA on grouped data
% 
% 
% 
% Tdisp = Tselect;
% 
% Tdisp = [Tdisp, table(fPeakLo), table(aPeakLo), table(avLoMag), ...
%     table(fPeakHi), table(aPeakHi), table(avHiMag), table(R), table(Pc)];
% 
% Tdisp = [Tdisp, table(accPreRawMean), table(accPreRawStdv), ...
%     table(accPrePercMean), table(accPrePercStdv), ...
%     table(accPreDbMean), table(accPreDbStdv), ...
%     table(accDbsRawMean), table(accDbsRawStdv), ...
%     table(accDbsPercMean), table(accDbsPercStdv), ...
%     table(accDbsDbMean), table(accDBSDbStdv), ...
%     table(accPosRawMean), table(accPosRawStdv), ...
%     table(accPosPercMean), table(accPosPercStdv), ...
%     table(accPosDbMean), table(accPosDbStdv)]; 
%     





% %% Display Peak frequency over time
% 
% tSinceHar = minutes(Tdisp.harRefTime(:)); % minutes
% 
% figure;
% % For the low band:
% subplot(2,4,1); plot(tSinceHar, Tdisp.fPeakLo, '-o'); 
% xlabel('Time since injection (mins)'); ylabel('Frequency (Hz)'); 
% title('Low band peak frequency (PSD)'); grid on
% 
% subplot(2,4,2); stem(tSinceHar, Tdisp.aPeakLo);
% xlabel('Time since injection (mins)'); ylabel('Power, A.U.'); 
% title('Low band peak Power (avPSD)'); grid on
% ax = gca; ax.YLim(2) = 40;
% 
% subplot(2,4,3); stem(tSinceHar, Tdisp.avLoMag);
% xlabel('Time since injection (mins)'); ylabel('Power, A.U.'); 
% title('Low band peak Power (avHilbert)'); grid on
% ax = gca; ax.YLim(2) = 0.1;
% 
% subplot(2,4,4); stem(tSinceHar, Tdisp.R);
% xlabel('Time since injection (mins)'); ylabel('Pearsons R'); 
% title('low-high band correlation'); grid on
% 
% 
% % For the high band: 
% subplot(2,4,5); plot(tSinceHar, Tdisp.fPeakHi, '-o'); 
% xlabel('Time since injection (mins)'); ylabel('Frequency (Hz)'); 
% title('High band peak frequency (PSD)'); grid on
% 
% subplot(2,4,6); stem(tSinceHar, Tdisp.aPeakHi);
% xlabel('Time since injection (mins)'); ylabel('Power, A.U.'); 
% title('High band peak Power (avPSD)'); grid on
% ax = gca; ax.YLim(2) = 40;
% 
% subplot(2,4,7); stem(tSinceHar, Tdisp.avHiMag);
% xlabel('Time since injection (mins)'); ylabel('Power, A.U.'); 
% title('High band peak Power (avHilbert)'); grid on
% ax = gca; ax.YLim(2) = 0.1;
% 
% 
% 
