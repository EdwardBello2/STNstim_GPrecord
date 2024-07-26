% batch for generating lots of iterations of the washin analysis

% NOTE: this script will only work on harmaline-day 

% TO-DO
% - test script as it currently is
% - test ability to Cache
% - refactor

%% Call toolboxes

%% Specify paths needed for inputs and outputs

% global PROJROOTPATH
% Add (temporarily) folder & subfolders of Git version-controlled code for Harmaline project
addpath(genpath('C:\Users\bello043\Documents\GitHub\Thalamic-DBS-HarmalineTremor'));

% Local path to google drive project folder
PROJROOTPATH = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';

% % Local NWB output 
% NWBPATH = [PROJROOTPATH 'Data Processing\NWB\'];



%% DEFINE all combinations of parameters for analysis as table with unique rows

dose = {'2 mg/kg', '4 mg/kg', '6 mg/kg', '8 mg/kg'};
expos = {'1st', '2nd'};
% winsec = 60; % seconds
lfpband = {[2, 4], [4, 8], [8, 12], [12, 35], [35, 55]}; % Hz
chanpair = {1, 2, 3, 4, 5, 6, 7};

B = allcomb(dose, expos, lfpband, chanpair);

dose = B(:,1);
expos = B(:,2);
lfpband = B(:,3);
chanpair = B(:,4);

Tp = [table(dose), table(expos), table(lfpband), table(chanpair)];



%% ITERATE the analysis for each row

WINSEC = 60; % seconds

nIter = height(Tp);
for i = 1:nIter
    close all
    func_harmalineWashin_NWB(cfg, Tp.dose{i}, Tp.expos{i}, WINSEC, Tp.lfpband{i}, ...
        Tp.chanpair{i}, PROJROOTPATH);
    disp('')
    
end


% 
% % Channel Pair strings
% for i = 1:8, cStr(i) = {['C' num2str(i-1)]}; end
% for i = 1:7, cpStr(i) = {[cStr{i} '-' cStr{i+1}]}; end
% 
% 
% 
% 
% cols.harmalineDose = DOSE;
% cols.exposure = EXPOS;
% cols.sessionType = 'exp'; % str, 'baseline' | 'exp'
% 
% T = readtable([NWBPATH 'nwbfile_epochIntervals_HHMMSS_exp.txt'], ...
%     'HeaderLines', 0);
% 
% subT = getRows(T, cols);
% nwbfn = subT.nwbFile_name{:}; % should only be one entry..
% 
% 
% 
% %% Initialize directories and filenames for intermediate (cache) data
% 
% % Make sure that intermediate data folder exists for this script
% intPath = assertIntermDir(PROJROOTPATH);
% 
% 
% % check if the data from a previous analysis has been cached for this
% % particular set of parameters
% namepars.nwbfn = nwbfn;
% namepars.cpStr = cpStr;
% namepars.LFP_BAND = LFP_BAND;
% namepars.intPath = intPath;
% namepars.cpair = cpair;
% namepars.WINSEC = WINSEC;
% 
% fullCachefn = genCacheName(namepars);
% isCached = checkCached(fullCachefn);
% 
% 
% 
% %% Load LFP data from nwb file
% 
% if ~isCached
%     % NOTE: use ignorecache option, believe me...
%     nwb = nwbRead([NWBPATH 'TremoLfpDBS-190927-100155.nwb'], 'ignorecache');
%     elecSeries = nwb.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('lfp');
% 
%     
%     % get timestamps
%     tst = elecSeries.timestamps.load;
%     lfp = elecSeries.data.load;
%     fs = elecSeries.starting_time_rate;
% 
% end
% 
% 
% 
% %% Filter for band of interest
% 
% if ~isCached
%     fclfp = LFP_BAND;
%     [blfp, alfp] = butter(2, fclfp / (fs/2), 'bandpass');
% 
%     lfpFilt = filtfilt(blfp, alfp, double(lfp));
%     
% end
% 
% 
% 
% %% Get subselection of LFP data as conntinuous timestamps and lfp data
% 
% if ~isCached
%     % Re-reference clock times of epochs to session_strattime
%     Tref = subT;
%     refCol = 5; % column for sessionStartType duration object
%     Tref{:,5:end} = subT{:,5:end} - subT{:,refCol};
% 
%     harDeliv = seconds(Tref.harDelivered(:));
%     harDeliv = harDeliv/60; % convert to minutes
% 
%     % pre-harmaline naive epoch
%     t_naivePre = [seconds(Tref.session_startTime(:)), seconds(Tref.naivePre_end(:))];
% 
%     % washin epoch
%     t_washin = [seconds(Tref.harDelivered(:)), seconds(Tref.washin_end(:))];
% 
%     % washout epoch
%     t_washout = [seconds(Tref.washout_beg(:)), seconds(Tref.session_endTime(:))];
% 
% 
%     % Get indices in data
%     isNaivePre = (tst >= t_naivePre(1)) & (tst < t_naivePre(2));
%       isWashin = (tst >= t_washin(1))   & (tst < t_washin(2));
%      isWashout = (tst >= t_washout(1))  & (tst < t_washout(2));
% 
%     isWashinExp = isNaivePre | isWashin | isWashout;
% 
%     
%     % Get sub-selection of LFP data and timestamps
%     lfp_naivePre = lfpFilt(isNaivePre,cpair);
%     tst_naivePre = tst(isNaivePre);
% 
%     lfp_washin = lfpFilt(isWashin,cpair);
%     tst_washin = tst(isWashin);
% 
%     lfp_washout = lfpFilt(isWashout,cpair);
%     tst_washout = tst(isWashout);
%     
% end
% 
% % figure; plot(tstSel, lfpSel(:,1));
%  
% 
% 
% %% Bin Times and RMS data
% 
% if ~isCached
%     % first get subselection of timeseries so that all bins ahve exactly 60
%     % seconds of data in them
%     winSamp = floor(WINSEC * fs);
% 
%     % Plot RMS over time, in 60 sec segments of signal
% 
%     disp('BEGIN RMS SLIDING FILTER...')
%     tic
%     rms_naivePre = smoothrms(lfp_naivePre, winSamp);
%     rms_washin = smoothrms(lfp_washin, winSamp);
%     rms_washout = smoothrms(lfp_washout, winSamp);
% 
%     rms_WashExp = [rms_naivePre; rms_washin; rms_washout];
% %     rms_WashExpnorm = rms_WashExp / mean(rms_naivePre);
%     rms_AvBaseline = mean(rms_naivePre);
% 
%     tst_WashExp = [tst_naivePre; tst_washin; tst_washout];
%     tst_WashExpmin = tst_WashExp / 60; % convert to mins
%     
%     tstmin = tst / 60;
%     disp('DONE!!!')
%     toc
% 
% end
% 
% 
% 
% %% Save analysis results into cache if not done yet, 
% % or load the cached data if it was already analyzed before
% 
% if ~isCached
%     save(fullCachefn, 'tstmin', 'harDeliv', 'tst_WashExpmin', 'rms_WashExp', 'rms_AvBaseline')
%      
% else
%     load(fullCachefn);
% 
% end
% 
% 
% %% PLOT lfp band over time
% 
% rms_WashExpnorm = rms_WashExp / rms_AvBaseline;
% 
% figure; ax = axes;
% % idx = 1:length(isWashinExp);
% % isIdx = idx(isWashinExp);
% % plot(tst_All, rms_All, 'MarkerIndices', isIdx);
% plot(tst_WashExpmin, rms_WashExpnorm);
% hold on
% plot([harDeliv, harDeliv], [ax.YLim(1), ax.YLim(2)], 'r')
% plot([max(tst_WashExpmin), max(tst_WashExpmin)], [ax.YLim(1), ax.YLim(2)], '--r')
% 
% ylabel('Norm to nairvePre, A.U.');
% xlabel('Time (min)');
% 
% titStr = [cpStr{cpair} ' ' EXPOS ' ' DOSE ' peri-harmaline [' num2str(LFP_BAND(1)) ' - ' ...
%     num2str(LFP_BAND(2)) ' Hz] band-power'];
% title(titStr);
% 
% 
% % 
% % % Get RMS values for each bin
% % nBins = floor(numel(washBand) / winSamp);
% % totSamps = nBins * winSamp;
% % washCorr = washBand(1:totSamps);
% % washBinned = reshape(washCorr, [winSamp, nBins]);
% % washRMS = rms(washBinned);
% % tBinsWash = winSec * (1:nBins);
% % 
% % 
% % % Get RMS values for baseline too
% % nBinsBase = floor(numel(baseBand) / winSamp);
% % totSampsBase = nBinsBase * winSamp;
% % baseCorr = baseBand(1:totSampsBase);
% % baseBinned = reshape(baseCorr, [winSamp, nBinsBase]);
% % baseRMS = rms(baseBinned);
% % tBinsBase = winSec * (1:nBinsBase);
% % tBinsBase = fliplr(tBinsBase) * (-1);
% % 
% % 
% % % Normalize by average baseline
% % baseRMSav = mean(baseRMS);
% % baseRMSnorm = baseRMS / baseRMSav;
% % washRMSnorm = washRMS / baseRMSav;
% % 
% % 
% % % Plot them together in time
% % titStr = [DAYREP,' ', DOSE, ' peri-harmaline [', num2str(fcband(1)), ' - ', ...
% %     num2str(fcband(2)), ' Hz] band-power'];
% % 
% % createWashinFig([tBinsBase, tBinsWash], [baseRMSnorm, washRMSnorm], titStr);
% 
% 
% 
% %% OLD CODE
% 
% % first-case: C0-C1
% 
% % % Recording of interest:
% % DATALOCPATH = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata';
% % 
% % % I pre-saved the constants below for each exp day, because I don't want to
% % % manually enter these things anew every frikkin time I want to check a
% % % frequency band power over time... look inside cfgExpDay func for the
% % % constants for each day:
% % % allowed strings: '1st_2mgkg' | '1st_4mgkg' | '1st_6mgkg' | '1st_8mgkg'
% % %                  '2nd_2mgkg' | '2nd_4mgkg' | '2nd_6mgkg' | '2nd_8mgkg'
% % [TDTTANK, TDTBLOCK, DOSE, ...
% %     DAYREP, washinBegClk, washinEndClk] = cfgExpDay('2nd_8mgkg');
% % 
% % 
% % % LFP frequency-range:
% % fclfp = [0.5, 300]; % Hz
% % 
% % 
% % % LFP band-power of interest: 
% % fcband = [12, 35]; % Hz
% % 
% % 
% % storeField_raw = 'RAW8'; 
% % chPair = [1, 2];
% % 
% % 
% % % Display parameters:
% % ELECS = {'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'};
% % DB_RANGE = [-170, -100];
% % FIGPOS = [2005 87 1387 859];
% % 
% % 
% % 
% % fulltdtpath = [DATALOCPATH, '\', TDTTANK, '\', TDTBLOCK];
% % 
% % tdt = TDTbin2mat(fulltdtpath, 'TYPE', {'epocs', 'snips', 'scalars'});
% % 
% % dt = datetime(tdt.info.date);
% % tt = datetime(tdt.info.utcStartTime); 
% % tott = dt + timeofday(tt);
% % sess_utcStart = datestr(tott, 'yyyy-mm-ddTHH:MM:SSZ');
% % 
% % 
% % 
% % %% choose data-limits for washin period
% % 
% % fileBegClk = datetime(tdt.info.Start, 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
% % 
% % washinBegT = washinBegClk - fileBegClk;
% % tbeg = seconds(washinBegT);
% % washinEndT = washinEndClk - fileBegClk;
% % tend = seconds(washinEndT);
% % 
% % 
% % 
% % %% load difference between channels, filter, downsample for lfp and track band-power over time
% % 
% % 
% % 
% % 
% % 
% % 
% % % Load washin data, and baseline data
% % washRaw = TDTbin2mat_chanDiff(fulltdtpath, chPair, 'STORE', storeField_raw, ...
% %     'RANGES', [tbeg; tend]);
% % 
% % tRangeBase = [0; 600]; % seconds
% % baseRaw = TDTbin2mat_chanDiff(fulltdtpath, chPair, 'STORE', storeField_raw, ...
% %     'RANGES', tRangeBase);
% % 
% % fs = washRaw.streams.(storeField_raw).fs;
% % 
% % 
% % 
% % % lfp-filter & downsample each
% % % fclfp = [0.5, 300];
% % [blfp, alfp] = butter(2, fclfp / (fs/2), 'bandpass');
% % 
% % washFlt = filtfilt(blfp, alfp, double(washRaw.streams.(storeField_raw).data));
% % baseFlt = filtfilt(blfp, alfp, double(baseRaw.streams.(storeField_raw).data));
% % 
% % dwnFactor = floor(fs/ 1000);
% % washLfp = downsample(washFlt, dwnFactor);
% % baseLfp = downsample(baseFlt, dwnFactor);
% % fsLfp = fs / dwnFactor;
% % 
% % 
% % 
% % % filter line noise from each
% % washClean = filtLineNoise(washLfp, fsLfp, 3);
% % baseClean = filtLineNoise(baseLfp, fsLfp, 3);
% % 
% % 
% % 
% % %% track power in 4-8 hz band for each
% % 
% % % fcband = [4, 8]; % hz
% % [b,a] = butter(4, fcband / (fsLfp/2), 'bandpass');
% % washBand = filtfilt(b, a, washClean);
% % baseBand = filtfilt(b, a, baseClean);
% % tLfp = (1/fsLfp) * (0:(length(washBand)-1));
% % 
% % 
% % % Plot RMS over time, in 60 sec segments of signal
% % 
% % % first get subselection of timeseries so that all bins ahve exactly 60
% % % seconds of data in them
% % binSec = 60; % seconds
% % binSamp = floor(binSec * fsLfp);
% % 
% % 
% % % Get RMS values for each bin
% % nBins = floor(numel(washBand) / binSamp);
% % totSamps = nBins * binSamp;
% % washCorr = washBand(1:totSamps);
% % washBinned = reshape(washCorr, [binSamp, nBins]);
% % washRMS = rms(washBinned);
% % tBinsWash = binSec * (1:nBins);
% % 
% % 
% % % Get RMS values for baseline too
% % nBinsBase = floor(numel(baseBand) / binSamp);
% % totSampsBase = nBinsBase * binSamp;
% % baseCorr = baseBand(1:totSampsBase);
% % baseBinned = reshape(baseCorr, [binSamp, nBinsBase]);
% % baseRMS = rms(baseBinned);
% % tBinsBase = binSec * (1:nBinsBase);
% % tBinsBase = fliplr(tBinsBase) * (-1);
% % 
% % 
% % % Normalize by average baseline
% % baseRMSav = mean(baseRMS);
% % baseRMSnorm = baseRMS / baseRMSav;
% % washRMSnorm = washRMS / baseRMSav;
% % 
% % 
% % % Plot them together in time
% % titStr = [DAYREP,' ', DOSE, ' peri-harmaline [', num2str(fcband(1)), ' - ', ...
% %     num2str(fcband(2)), ' Hz] band-power'];
% % 
% % createWashinFig([tBinsBase, tBinsWash], [baseRMSnorm, washRMSnorm], titStr);
% % 
% % disp('')
% 
% 
% 
% %% SUB-FUNCTIONS
% 
% function [isCached] = checkCached(fullPathFn)
% isCached = false;
% 
% 
% if exist(fullPathFn, 'file')
%     isCached = true;
% %     load(fullPathFn)
% end
% 
% end
% 
% function fullPathFn = genCacheName(namepars)
% 
% nwbfn = namepars.nwbfn;
% cpStr = namepars.cpStr;
% LFP_BAND = namepars.LFP_BAND;
% intPath = namepars.intPath;
% cpair = namepars.cpair;
% WINSEC = namepars.WINSEC;
% 
% [~,nwblab,~] = fileparts(nwbfn);
% 
% formatSpecFn = 'RMS_%s_AvWin%ss_ch%s_lfpHz%s';
% matfnStr = sprintf(formatSpecFn, ...
%     nwblab, ...
%     num2str(WINSEC), ...
%     cpStr{cpair}, ...
%     [num2str(LFP_BAND(1)) '-' num2str(LFP_BAND(2))]);
% 
% fullPathFn = [intPath matfnStr '.mat'];
% 
% end
% 
% function interDataPath = assertIntermDir(PROJROOTPATH)
% % Get the name of the currently running script
% 
% % global PROJROOTPATH
% 
% [scriptDirectoryFullPath, scriptName] = fileparts(mfilename('fullpath'));
% 
% % First make sure that this script has a folder within the project folder's
% % intermediate data section
% interDataDir = [PROJROOTPATH 'Data Processing\intermediateData\' scriptName];
% if ~exist(interDataDir, 'dir')
%     mkdir(interDataDir)
%     
% end
% 
% interDataPath = [interDataDir '\'];
% 
% end
% 
% function xrms = smoothrms(x, winSamps)
% % First get the square of every point, then smooth with a moving average
% % square window, finally get the square root of every point.
% 
% w = ones(winSamps,1);
% w = w ./ sum(w);
% 
% xSqr = x.^2;
% xSqrSmooth = conv(xSqr, w, 'same');
% 
% xrms = sqrt(xSqrSmooth);
% 
% end
% 
% function subT = getRows(T, cols)
% 
% colNames = fieldnames(cols);
% nFields = length(fieldnames(cols));
% 
% isIdx = false(height(T), nFields);
% for iF = 1:nFields
%     isIdx(:,iF) = strcmp(cols.(colNames{iF}), T{:,colNames{iF}});
%     
% end
% 
% select = all(isIdx, 2);
% subT = T(select,:);
% 
% 
% end
% 
% function createWashinFig(X1, Y1, titStr)
% %CREATEFIGURE(X1, Y1, X2, Y2)
% %  X1:  vector of x data
% %  Y1:  vector of y data
% %  X2:  vector of x data
% %  Y2:  vector of y data
% 
% %  Auto-generated by MATLAB on 27-Feb-2020 11:46:14
% 
% % Create figure
% figure1 = figure;
% 
% % Create axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% 
% % Create plot
% plot(X1,Y1,'Marker','o','LineStyle','none');
% 
% % Create plot
% plot([0, 0], [axes1.YLim(1), axes1.YLim(2)], '--');
% 
% % Create ylabel
% ylabel('Power norm to baseline (A.U.)');
% 
% % Create xlabel
% xlabel('Time (s)');
% 
% % Create title
% title(titStr);
% 
% box(axes1,'on');
% grid(axes1,'on');
% 
% end
% 
% function [ax, fig] = getmtspecgramfig(chans, blockPath)
% % function generates a figure of the multitaper spectrogram resulting from
% % a bipolar lfp. Bipolar lfp comes by extracting two LFP data streams from
% % the TDT file according to the channel-pair specified in "chans".
% 
% 
% % extract 1 channel and filter/down-sample it
% tic
% % extract1 = TDTbin2mat(blockPath, 'RANGES', [s1; s2], 'STORE', 'RAW8', 'CHANNEL', chans(1));
% extract1 = TDTbin2mat(blockPath, 'STORE', 'RAW8', 'CHANNEL', chans(1));
% sig1 = extract1.streams.RAW8.data';
% fs = extract1.streams.RAW8.fs;
% toc
% disp('done loading first chan')
% clear extract1
% 
% 
% tic
% % extract2 = TDTbin2mat(blockPath, 'RANGES', [s1; s2], 'STORE', 'RAW8', 'CHANNEL', chans(2));
% extract2 = TDTbin2mat(blockPath, 'STORE', 'RAW8', 'CHANNEL', chans(2));
% sig2 = extract2.streams.RAW8.data';
% toc
% disp('done loading second chan')
% clear extract2
% 
% lfp = sig1 - sig2;
% clear sig1 sig2
% 
% tic
% % [lfpDwnSmp, fsDwnSmp] = preconditionSig(double(lfp), fs);
% fsDwnSmp = fs;
% lfpDwnSmp = double(lfp);
% toc
% 
% disp('done filtering/downsampling LFP')
% 
% 
% 
% %% Get re-referenced LFPs
% 
% % redo this section with Chronux methods instead
% 
% params.Fs = fs;
% fc = 100;
% % params.fpass = [0, 300];
% params.tapers = [3, 5];
% movingwin = [5, 1];
% 
% 
% % [S, f] = mtspectrumc(lfp, params);
% disp('begin spectrogram')
% tic
% [S, t, f] = mtspecgramc(lfp, movingwin, params);
% clear lfp
% % S = S';
% toc
% disp('end spectrogram')
% 
% % display only frequencies below 300Hz
% fIdx = f <= 300;
% fLim = f(fIdx);
% 
% % also getting spectrogram fft's to display on surface plot
% % sNorm = abs(s);
% fig = figure; ax = axes;
% surf(t, f(fIdx), 10*log10(S(:,fIdx)'+eps), 'EdgeColor', 'none');  
% axis xy; 
% axis tight; 
% colormap(parula); 
% view(0,90);
% % ax.YLim = [0, 100];
% % cBar = colorbar;
% % cBar.Label.String = 'Power/frequency (dB/Hz)';
% % caxis([-170, -100])
% % 
% % xlabel('Time (s)')
% % ylabel('Frequency (Hz)')
% % 
% % title([ELECS{chans(1)}, '-', ELECS{chans(2)}]);
% % fig.Position = [2005 87 1387 859];
% clear t f S
% 
% 
% end
% 
% function [TDTTANK,TDTBLOCK,DOSE,DAYREP,washinBegClk,washinEndClk] = cfgExpDay(expStr)
% 
% switch expStr
%     case '2nd_8mgkg'
%         TDTTANK = 'Uva-191115';
%         TDTBLOCK = 'TremoLfpDBS-191115-100127';
%         % YYYYMMDD = '20191018';
%         DOSE = '8mgkg';
%         DAYREP = '2nd';
% 
%         washinBegClk = datetime('10:32:50 2019-11-15', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         washinEndClk = datetime('11:33:04 2019-11-15', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         
%     case '2nd_6mgkg'
%         TDTTANK = 'Uva-191108';
%         TDTBLOCK = 'TremoLfpDBS-191108-101829';
%         % YYYYMMDD = '20191018';
%         DOSE = '6mgkg';
%         DAYREP = '2nd';
% 
%         washinBegClk = datetime('10:50:09 2019-11-08', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         washinEndClk = datetime('11:50:46 2019-11-08', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         
%     case '2nd_4mgkg'
%          TDTTANK = 'Uva-191101';
%         TDTBLOCK = 'TremoLfpDBS-191101-101430';
%         % YYYYMMDD = '20191018';
%         DOSE = '4mgkg';
%         DAYREP = '2nd';
% 
%         washinBegClk = datetime('10:48:17 2019-11-01', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         washinEndClk = datetime('11:48:34 2019-11-01', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         
%     case '2nd_2mgkg'
%         TDTTANK = 'Uva-191025';
%         TDTBLOCK = 'TremoLfpDBS-191025-104651';
%         % YYYYMMDD = '20191018';
%         DOSE = '2mgkg';
%         DAYREP = '2nd';
% 
%         washinBegClk = datetime('11:17:59 2019-10-25', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         washinEndClk = datetime('12:18:20 2019-10-25', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         
%     case '1st_8mgkg'
%         TDTTANK = 'Uva-191018';
%         TDTBLOCK = 'TremoLfpDBS-191018-100615';
%         % YYYYMMDD = '20191018';
%         DOSE = '8mgkg';
%         DAYREP = '1st';
% 
%         washinBegClk = datetime('10:38:45 2019-10-18', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         washinEndClk = datetime('11:39:37 2019-10-18', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         
%     case '1st_6mgkg'
%          TDTTANK = 'Uva-191011';
%         TDTBLOCK = 'TremoLfpDBS-191011-104322';
%         % YYYYMMDD = '20191018';
%         DOSE = '6mgkg';
%         DAYREP = '1st';
% 
%         washinBegClk = datetime('11:14:54 2019-10-11', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         washinEndClk = datetime('12:14:48 2019-10-11', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         
%     case '1st_4mgkg'
%         TDTTANK = 'Uva-191004';
%         TDTBLOCK = 'TremoLfpDBS-191004-100637';
%         % YYYYMMDD = '20191018';
%         DOSE = '4mgkg';
%         DAYREP = '1st';
% 
%         washinBegClk = datetime('10:39:59 2019-10-04', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         washinEndClk = datetime('11:40:03 2019-10-04', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         
%     case '1st_2mgkg'
%         TDTTANK = 'Uva-190927';
%         TDTBLOCK = 'TremoLfpDBS-190927-100155';
%         % YYYYMMDD = '20191018';
%         DOSE = '2mgkg';
%         DAYREP = '1st';
% 
%         washinBegClk = datetime('10:31:33 2019-09-27', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         washinEndClk = datetime('11:31:43 2019-09-27', 'InputFormat', 'HH:mm:ss yyyy-MM-dd');
%         
%     otherwise
%         error('wrong string entry! check cfgExpDay func at bottom of this script')
%         
% end
% 
% 
% end
