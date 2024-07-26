% script to read in and process Kiwi's tremor data, based on:
% "script_melaTremor_PSDs.m"

% Collected on microstrain G-link hardware
clear; close all

% G-link sampling freq when streaming 3 channels (hardware design)
% fs = 617; % samples/second



%% Load in data and get into format of "acc" and "fs"


% for loading
dataAcqPn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\';
directoryPn = 'Kiwi\KiwiExperiment_05_17_2016\';



% for saving
projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';
SavePn = 'Data Processing\Kiwi_tremors\';

% recTimes = xlsread([projRootPath 'Data Acquisition\Kiwi Experiment\' 'RecordingsMetadata' '.xlsx'], ...
%     'B2:B34');
% a = datetime(recTimes, 'convertfrom', 'excel');
% a.Format = 'HH:mm:ss'

% T = readtable([projRootPath 'Data Acquisition\Kiwi Experiment\' 'RecordingsMetadata2' '.xlsx'])
% a = datetime(T.recBeginText)
% har = datetime('12:13:00', 'InputFormat', 'HH:mm:ss')% Read in all data for .mat file

% data filename
fn = 'tremor_frequency_measure_05-17-16_time1220_arm';

load([dataAcqPn directoryPn fn '.mat'], 'data', 'time', 'sampling_frequency', 'output_signal')


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
    


%% visualize code

% View raw trace of accel data
figure; plot(t, data)
xlabel('Time (seconds)');

% detrend the data with highpass filter
fc = 1; 
[b,a] = butter(2, fc/(fs/2), 'high');
accFilt = filtfilt(b, a, acc);

figure; plot(t, accFilt)
title('highpassed')


% show spectrogram of data

% run spectrogram on it
tWin = 1; % seconds
window = floor(tWin * fs);
% noverlap = floor(0.5*window);
noverlap = 0;



% Sum PSDs from all 3 accel channels
for ch = 1:3
%     ch = 1;
    [Z(:,:,ch), f, tSpec] = spectrogram(accFilt(:,ch), window, noverlap, [], fs, 'yaxis');
    % tRef = linspace(tSpec(1), tSpec(end), numel(tSpec));
    S(:,:,ch) = abs(Z(:,:,ch));
    S(:,:,ch) = S(:,:,ch).^2;

end
% Combine spectral power from all three axes (summate)
Sacc = sum(S, 3);

% Get average power (S) PSD over time
% PSDWashNorm = Swash ./ (repmat(PSDbaseAv, 1, size(Swash, 2)));
sm = 500;
f2 = figure; 
ax = axes;
ax.Parent = f2;
surf(tSpec, f, 10*log10(Sacc), ...
    'Parent', ax, ...
    'EdgeColor', 'none');
hold(ax, 'on');
% xlabel('Time (min)');
grid on
axis(ax, 'tight')
view(0,90)
colorbar
% ax.CLim = [-9, 9];
% ax.YLim = [0, 20];
% title([titStr, ', Washin, dBgain (Washin/Naive)'])


%% Clean spectrogram of movement artifacts (sudden jerks)

[~, isMovArt] = util.remoutliers(sum(Sacc, 1), ...
    'bound', 'upper', ...
    'MADthresh', 2);
temp = Sacc;
temp(:,isMovArt) = [];
SaccMarkOutliers = Sacc;
SaccMarkOutliers(:,isMovArt) = 100;


% Show Spectrogram with detected movement artifacts marked as 100's values
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
title(['Mov Artifacts marked'])
ylabel('Frequency (Hz)');
f4.Position = [2555 57 560 421];


% Show Spectrogram with movement artifacts removed totally
f4 = figure;
ax4 = axes;
ax4.Parent = f4;
tset = t;
tset(isMovArt) = [];
surf([1:size(temp, 2)], f, 10*log10(temp), ...
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
title(['Mov Artifacts removed'])
ylabel('Frequency (Hz)');
f4.Position = [2554 566 560 420];



%% Display average "clean" PSD for the entire recording period

Sclean = temp; % movement artifact windows removed
PSD = mean(Sclean, 2);

f1 = figure; ax = axes;
plot(f, (PSD));
grid on; 
ax.XLim = [0, 35]
xlabel('Frequency (Hz)');
title(fn, 'Interpreter', 'none')



%% Fill out metadata fields into a table row

% time since harmaline injection


% sampling frequency
samplingFrequency = sampling_frequency;

% total time in recording, seconds
totTime = time(end);

% percentage of that time occupied by movement artifacts (1-sec window
% detection)
movArtPerc = sum(isMovArt) * tWin / totTime;

% max value in PSD, for help in scaling results later



%% Save results to Processing folder
[~,savefn,~] = fileparts(fn);


% save([projRootPath SavePn savefn '_PSD'], 'Sacc', 'tSpec', 'f', 'isMovArt', 'Sclean');
% savefig(f1, [projRootPath SavePn savefn '_PSD']);



%% % test cross-correlation between two signals
% 
% % figure; plot(acc(:,1))
% 
% % get low band and high band
% 
% fcLow = [2, 10]; % Hz
% fcHigh = [12, 20]; % Hz
% order = 2; 
% 
% [b,a] = butter(order, fcLow / (fs/2), 'bandpass');
% acc1Low = filtfilt(b, a, acc(:,1));
% 
% [b,a] = butter(order, fcHigh / (fs/2), 'bandpass');
% acc1High = filtfilt(b, a, acc(:,1));
% 
% figure; hold on
% plot(acc1Low); plot(acc1High);
% 
% 
% [r, lags] = xcorr(acc1Low, acc1High);
% 
% figure; plot(lags/fs, r); grid on
% 
% 



