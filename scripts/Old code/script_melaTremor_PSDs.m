% sample script to read in Mela data

% Collected on microstrain G-link hardware
clear; close all

% G-link sampling freq when streaming 3 channels (hardware design)
fs = 617; % samples/second

dataAcqPn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\';
directoryPn = 'Mela\Harmaline\';
projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';
SavePn = 'Data Processing\Mela_tremors\';



% data filename
fn = '09-28-2015_1056_PostHarmaline.csv';



% Read in accel data from a given .csv file
T = readtable([dataAcqPn directoryPn fn]);




% View raw trace of accel data
t = 1/fs * (0:(height(T) - 1))';
acc = [T.Channel1, T.Channel2, T.Channel3];
figure; plot(t, acc)
figure; plot(t, T.AcquisitionAttempt)
xlabel('Time (seconds)');
ylabel('Sample');
title('sample num vs. time of occurrence')

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
% ax4.YLim = [0, 300];
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


% Save results to Processing folder
[~,savefn,~] = fileparts(fn);


save([projRootPath SavePn savefn '_PSD'], 'Sacc', 'tSpec', 'f', 'isMovArt', 'Sclean');
savefig(f1, [projRootPath SavePn savefn '_PSD']);
