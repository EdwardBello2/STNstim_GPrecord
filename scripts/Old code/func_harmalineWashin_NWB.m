% batch for generating lots of iterations of the washin analysis

% NOTE: this script will only work on harmaline-day 

% TO-DO
% - test script as it currently is
% - test ability to Cache
% - refactor

%% Call toolboxes

function func_harmalineWashin_NWB(nwb, cfg, rmswin, lfpband, cpair, sessionMetaData)
%% Specify paths needed for inputs and outputs

% global PROJROOTPATH
% Add (temporarily) folder & subfolders of Git version-controlled code for Harmaline project
% addpath(genpath('C:\Users\bello043\Documents\GitHub\Thalamic-DBS-HarmalineTremor'));

% Local path to google drive project folder
% PROJROOTPATH = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';
% 
% % Local path to Original TDT data
% TDTLOCPATH = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\'; 
% 
% % Local path to LFP data
% LFPLOCPATH = [PROJROOTPATH 'Data Processing\LFPdownsamp\'];
% SUFFIX = '_diffLFP'; % append to TDTblock filename to get .mat file containing converted LFP data

% Local NWB output 
% NWBPATH = [cfg.projRootDirFullpath 'Data Processing\NWB\'];
DOSE = sessionMetaData.dose;
EXPOS = sessionMetaData.exposureNum;
SESSTYPE = sessionMetaData.sessionType;

% isCached = false; % default for tracking whether or not a given run of the script has been done before


%% Read in data of interest from NWB

% global DOSE EXPOS LFP_BAND cpair cStr cpStr nwbfn intPath WINSEC

% get subselection from table for a given session
% DOSE = '2 mg/kg'; % str, '2 mg/kg' | '4 mg/kg' | '6 mg/kg' | '8 mg/kg'
% EXPOS = '1st'; % str, '1st' | '2nd'
% WINSEC = 60; % seconds
% LFP_BAND = [4, 8]; % Hz

% cpair = 1;

% Channel Pair strings
for i = 1:8, cStr(i) = {['C' num2str(i-1)]}; end
for i = 1:7, cpStr(i) = {[cStr{i} '-' cStr{i+1}]}; end




% cols.harmalineDose = DOSE;
% cols.exposure = EXPOS;
% % cols.sessionType = 'exp'; % str, 'baseline' | 'exp'
% cols.sessionType = SESSTYPE;
% 
% T = readtable([NWBPATH 'nwbfile_epochIntervals_HHMMSS_exp.txt'], ...
%     'HeaderLines', 0);
% 
% subT = getRows(T, cols);
% nwbfn = subT.nwbFile_name{:}; % should only be one entry..
% 
% [~, nwbfn, ~] = fileparts(nwb.general_session_id);
% 
% %% Initialize directories and filenames for intermediate (cache) data
% 
% % Make sure that cache data folder exists for this script
% [~, scriptName, ~] = fileparts(mfilename('fullpath'));
% cacheFullpath = assertCacheDir(cfg, scriptName);
% 
% 
% % check if the data from a previous analysis has been cached for this
% % particular set of parameters
% namepars.nwbfn = nwbfn;
% namepars.cpStr = cpStr;
% namepars.LFP_BAND = LFP_BAND;
% namepars.intPath = cacheFullpath;
% namepars.cpair = cpair;
% namepars.WINSEC = WINSEC;
% 
% fullCachefn = genCacheName(namepars);
% isCached = checkCached(fullCachefn);
% 
% 
% 
% % Load LFP data from nwb file
% if ~isfield(sessionMetaData, 'epochs')
%     sessionMetaData.epochs = {'all'};
% 
% end

% For each epoch-string specified in sessionMetaData.epochs:
nEps = numel(sessionMetaData.epochs);
for iEp = 1:nEps
    
    [lfprmsCell{iEp},tstCell{iEp},fs] = calcrms_nwblfp(nwb, rmswin, cfg, ...
        'epoch', sessionMetaData.epochs{iEp}, ...
        'channel', cpair, ...
        'lfpband', lfpband);
    
    
    
%     % Get lfps
%     [lfp, tst, fs] = nwbLoadLfp(nwb, ...
%         'channel', 1, ...
%         'epoch', sessionMetaData.epochs{iEp});
% 
% 
%     % Filter for band of interest
%     fclfp = LFP_BAND;
%     [blfp, alfp] = butter(2, fclfp / (fs/2), 'bandpass');
% 
% %     lfpFilt = filtfilt(blfp, alfp, double(lfp));
%     lfpFilt = bandfilt_edgecomp(blfp, alfp, double(lfp));
%     
% 
%     % Bin Times and RMS data
%     % first get subselection of timeseries so that all bins ahve exactly 60
%     % seconds of data in them
%     winSamp = floor(WINSEC * fs);
% 
%     % Plot RMS over time, in 60 sec segments of signal
% 
%     disp('BEGIN RMS SLIDING FILTER...')
%     tic
% %     rms_washin = smoothrms(lfp_washin, winSamp);
% %     rms_washout = smoothrms(lfp_washout, winSamp);
% 
% %     rms_WashExp = [rms_naivePre; rms_washin; rms_washout];
% %     rms_AvBaseline = mean(rms_naivePre);
% 
% %     tst_WashExp = [tst_naivePre; tst_washin; tst_washout];
% %     tst_WashExpmin = tst_WashExp / 60; % convert to mins
% 
% %     tstmin = tst / 60;
% 
%     lfpRmsCell{iEp} = smoothrms(lfpFilt, winSamp);
%     tstminCell{iEp} = tst / 60;
%     disp('DONE!!!')
%     toc
    

end


% Combine individual epochs data into one vector/matrix
rms_WashExp = [];
tst = [];
for iEp = 1:nEps
    if strcmp('naivePreExp', sessionMetaData.epochs{iEp}) % if it's the baseline reference epoch
        % keep only the first 10 minutes of reference epoch (600 sec)
        rms_Baseline = lfprmsCell{iEp};
        tst_Baseline = tstCell{iEp};

        isfirstTen = tst_Baseline < 600;
        rms_Baseline = rms_Baseline(isfirstTen);
        tst_Baseline = tst_Baseline(isfirstTen);
        
        lfprmsCell{iEp} = rms_Baseline;
        tstCell{iEp} = tst_Baseline;
        
    end
        
    
    rms_WashExp = [rms_WashExp; lfprmsCell{iEp}];
    tst = [tst; tstCell{iEp}];
        
end

% Normalize rms by baseline condition
% get first 10 mins of data for normalization (large transients present
% when switching electrodes after 10 min
tstmin = tst / 60; 
% tstmin = tst; 

rms_AvBaseline = mean(lfprmsCell{1}); % initial naive condition
rms_WashExpnorm = rms_WashExp / rms_AvBaseline;



%% PLOT lfp band over time

% rms_WashExpnorm = rms_WashExp / rms_AvBaseline;

figure; ax = axes;
% idx = 1:length(isWashinExp);
% isIdx = idx(isWashinExp);
% plot(tst_All, rms_All, 'MarkerIndices', isIdx);
plot(tstmin, rms_WashExpnorm);
hold on
% plot([harDeliv, harDeliv], [ax.YLim(1), ax.YLim(2)], 'r')
% plot([max(tst_WashExpmin), max(tst_WashExpmin)], [ax.YLim(1), ax.YLim(2)], '--r')

ylabel('Norm to naivePre, A.U.');
xlabel('Time (min)');

titStr = [cpStr{cpair} ' ' EXPOS ' ' DOSE ' peri-harmaline [' num2str(lfpband(1)) ' - ' ...
    num2str(lfpband(2)) ' Hz] band-power'];
title(titStr);

end



%% SUB-FUNCTIONS

function xFilt = bandfilt_edgecomp(b, a, x)
% run filtfilt while compensating for edge effects by padding the end with
% flipped-mirrored end-portions prefixed and appended. x is col vector of
% data. 

R = 0.1; % 10% of signal
Nr = 50; % max num of points
N = size(x,1);
NR = min(round(N*R),Nr); % At most 50 points
for i = 1:size(x,2)
    x1(:,i) = 2*x(1,i) - flipud(x(2:NR+1,i));  % maintain continuity in level and slope
    x2(:,i) = 2*x(end,i) - flipud(x(end-NR:end-1,i));
end
x = [x1;x;x2];
% Do filtering
xFilt = filtfilt(b, a, x);
xFilt = xFilt(NR+1:end-NR,:);


end

function [isCached] = checkCached(fullPathFn)
isCached = false;


if exist(fullPathFn, 'file')
    isCached = true;
%     load(fullPathFn)
end

end

function fullPathFn = genCacheName(namepars)

nwbfn = namepars.nwbfn;
cpStr = namepars.cpStr;
LFP_BAND = namepars.LFP_BAND;
intPath = namepars.intPath;
cpair = namepars.cpair;
WINSEC = namepars.WINSEC;

[~,nwblab,~] = fileparts(nwbfn);

formatSpecFn = 'RMS_%s_AvWin%ss_ch%s_lfpHz%s';
matfnStr = sprintf(formatSpecFn, ...
    nwblab, ...
    num2str(WINSEC), ...
    cpStr{cpair}, ...
    [num2str(LFP_BAND(1)) '-' num2str(LFP_BAND(2))]);

fullPathFn = [intPath matfnStr '.mat'];

end

function interDataPath = assertIntermDir(PROJROOTPATH)
% Get the name of the currently running script

% global PROJROOTPATH

[scriptDirectoryFullPath, scriptName] = fileparts(mfilename('fullpath'));

% First make sure that this script has a folder within the project folder's
% intermediate data section
interDataDir = [PROJROOTPATH 'Data Processing\intermediateData\' scriptName];
if ~exist(interDataDir, 'dir')
    mkdir(interDataDir)
    
end

interDataPath = [interDataDir '\'];

end

function xrms = smoothrms(x, winSamps)
% First get the square of every point, then smooth with a moving average
% square window, finally get the square root of every point.

w = ones(winSamps,1);
w = w ./ sum(w);

xSqr = x.^2;
xSqrSmooth = conv(xSqr, w, 'same');

xrms = sqrt(xSqrSmooth);

end

function subT = getRows(T, cols)

colNames = fieldnames(cols);
nFields = length(fieldnames(cols));

isIdx = false(height(T), nFields);
for iF = 1:nFields
    isIdx(:,iF) = strcmp(cols.(colNames{iF}), T{:,colNames{iF}});
    
end

select = all(isIdx, 2);
subT = T(select,:);


end

function createWashinFig(X1, Y1, titStr)
%CREATEFIGURE(X1, Y1, X2, Y2)
%  X1:  vector of x data
%  Y1:  vector of y data
%  X2:  vector of x data
%  Y2:  vector of y data

%  Auto-generated by MATLAB on 27-Feb-2020 11:46:14

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
plot(X1,Y1,'Marker','o','LineStyle','none');

% Create plot
plot([0, 0], [axes1.YLim(1), axes1.YLim(2)], '--');

% Create ylabel
ylabel('Power norm to baseline (A.U.)');

% Create xlabel
xlabel('Time (s)');

% Create title
title(titStr);

box(axes1,'on');
grid(axes1,'on');

end

function [ax, fig] = getmtspecgramfig(chans, blockPath)
% function generates a figure of the multitaper spectrogram resulting from
% a bipolar lfp. Bipolar lfp comes by extracting two LFP data streams from
% the TDT file according to the channel-pair specified in "chans".


% extract 1 channel and filter/down-sample it
tic
% extract1 = TDTbin2mat(blockPath, 'RANGES', [s1; s2], 'STORE', 'RAW8', 'CHANNEL', chans(1));
extract1 = TDTbin2mat(blockPath, 'STORE', 'RAW8', 'CHANNEL', chans(1));
sig1 = extract1.streams.RAW8.data';
fs = extract1.streams.RAW8.fs;
toc
disp('done loading first chan')
clear extract1


tic
% extract2 = TDTbin2mat(blockPath, 'RANGES', [s1; s2], 'STORE', 'RAW8', 'CHANNEL', chans(2));
extract2 = TDTbin2mat(blockPath, 'STORE', 'RAW8', 'CHANNEL', chans(2));
sig2 = extract2.streams.RAW8.data';
toc
disp('done loading second chan')
clear extract2

lfp = sig1 - sig2;
clear sig1 sig2

tic
% [lfpDwnSmp, fsDwnSmp] = preconditionSig(double(lfp), fs);
fsDwnSmp = fs;
lfpDwnSmp = double(lfp);
toc

disp('done filtering/downsampling LFP')



%% Get re-referenced LFPs

% redo this section with Chronux methods instead

params.Fs = fs;
fc = 100;
% params.fpass = [0, 300];
params.tapers = [3, 5];
movingwin = [5, 1];


% [S, f] = mtspectrumc(lfp, params);
disp('begin spectrogram')
tic
[S, t, f] = mtspecgramc(lfp, movingwin, params);
clear lfp
% S = S';
toc
disp('end spectrogram')

% display only frequencies below 300Hz
fIdx = f <= 300;
fLim = f(fIdx);

% also getting spectrogram fft's to display on surface plot
% sNorm = abs(s);
fig = figure; ax = axes;
surf(t, f(fIdx), 10*log10(S(:,fIdx)'+eps), 'EdgeColor', 'none');  
axis xy; 
axis tight; 
colormap(parula); 
view(0,90);
% ax.YLim = [0, 100];
% cBar = colorbar;
% cBar.Label.String = 'Power/frequency (dB/Hz)';
% caxis([-170, -100])
% 
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% 
% title([ELECS{chans(1)}, '-', ELECS{chans(2)}]);
% fig.Position = [2005 87 1387 859];
clear t f S


end

function [TDTTANK,TDTBLOCK,DOSE,DAYREP,washinBegClk,washinEndClk] = cfgExpDay(expStr)

switch expStr
    case '2nd_8mgkg'
        TDTTANK = 'Uva-191115';
        TDTBLOCK = 'TremoLfpDBS-191115-100127';
        % YYYYMMDD = '20191018';
        DOSE = '8mgkg';
        DAYREP = '2nd';

        washinBegClk = datetime('10:32:50 2019-11-15', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        washinEndClk = datetime('11:33:04 2019-11-15', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        
    case '2nd_6mgkg'
        TDTTANK = 'Uva-191108';
        TDTBLOCK = 'TremoLfpDBS-191108-101829';
        % YYYYMMDD = '20191018';
        DOSE = '6mgkg';
        DAYREP = '2nd';

        washinBegClk = datetime('10:50:09 2019-11-08', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        washinEndClk = datetime('11:50:46 2019-11-08', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        
    case '2nd_4mgkg'
         TDTTANK = 'Uva-191101';
        TDTBLOCK = 'TremoLfpDBS-191101-101430';
        % YYYYMMDD = '20191018';
        DOSE = '4mgkg';
        DAYREP = '2nd';

        washinBegClk = datetime('10:48:17 2019-11-01', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        washinEndClk = datetime('11:48:34 2019-11-01', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        
    case '2nd_2mgkg'
        TDTTANK = 'Uva-191025';
        TDTBLOCK = 'TremoLfpDBS-191025-104651';
        % YYYYMMDD = '20191018';
        DOSE = '2mgkg';
        DAYREP = '2nd';

        washinBegClk = datetime('11:17:59 2019-10-25', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        washinEndClk = datetime('12:18:20 2019-10-25', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        
    case '1st_8mgkg'
        TDTTANK = 'Uva-191018';
        TDTBLOCK = 'TremoLfpDBS-191018-100615';
        % YYYYMMDD = '20191018';
        DOSE = '8mgkg';
        DAYREP = '1st';

        washinBegClk = datetime('10:38:45 2019-10-18', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        washinEndClk = datetime('11:39:37 2019-10-18', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        
    case '1st_6mgkg'
         TDTTANK = 'Uva-191011';
        TDTBLOCK = 'TremoLfpDBS-191011-104322';
        % YYYYMMDD = '20191018';
        DOSE = '6mgkg';
        DAYREP = '1st';

        washinBegClk = datetime('11:14:54 2019-10-11', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        washinEndClk = datetime('12:14:48 2019-10-11', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        
    case '1st_4mgkg'
        TDTTANK = 'Uva-191004';
        TDTBLOCK = 'TremoLfpDBS-191004-100637';
        % YYYYMMDD = '20191018';
        DOSE = '4mgkg';
        DAYREP = '1st';

        washinBegClk = datetime('10:39:59 2019-10-04', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        washinEndClk = datetime('11:40:03 2019-10-04', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        
    case '1st_2mgkg'
        TDTTANK = 'Uva-190927';
        TDTBLOCK = 'TremoLfpDBS-190927-100155';
        % YYYYMMDD = '20191018';
        DOSE = '2mgkg';
        DAYREP = '1st';

        washinBegClk = datetime('10:31:33 2019-09-27', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        washinEndClk = datetime('11:31:43 2019-09-27', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
        
    otherwise
        error('wrong string entry! check cfgExpDay func at bottom of this script')
        
end


end
