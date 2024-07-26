% Reads in nwb file and user parameters specifying lfp data to convert to
% rms signal power

function [lfpRMS, timestamps, fs] = calcrms_nwblfp(nwb, varargin)

% calcrms_nwblfp(nwb, rmswin, epoch, channel, lfpband, cfg)

% INPUTS
% DOSE = sessionMetaData.dose;
% EXPOS = sessionMetaData.exposureNum;

defaultEpoch = 'all';
defaultChan = [];
defaultBand = [];
defaultCfg = [];

validScalar = @(x) isnumeric(x) && isscalar(x);
% validIdx = @(x) size(x, 1) isisinteger(int8(x))
% validVector = @(x) isinteger(x) && isvector(x);
validBand = @(x) isnumeric(x) && isvector(x) && (length(x) < 3);

p = inputParser;
 addRequired(p, 'nwb');
 addRequired(p, 'rmswin', validScalar);
addParameter(p, 'epoch', defaultEpoch, @ischar);
addParameter(p, 'channel', defaultChan, @isvector);
addParameter(p, 'lfpband', defaultBand, validBand);
 addOptional(p, 'cfg', defaultCfg, @isstruct);
parse(p, nwb, varargin{:});


     nwb = p.Results.nwb;
  rmswin = p.Results.rmswin;
      ch = p.Results.channel;
epochStr = p.Results.epoch;
 lfpband = p.Results.lfpband;
     cfg = p.Results.cfg;
     

     
%% Caching functionality

if isempty(cfg)
    cfg.allowCache = false; 
    cfg.overwriteCache = false;
    
end

if cfg.allowCache
    % Generate unique name for data results in this run
    [~, nwbName, ~] = fileparts(nwb.general_session_id);
    matfnStr = [nwbName, 'rmswin', num2str(rmswin), epochStr, ...
        'ch', sprintf('-%d', ch), 'lfp', sprintf('-%d', lfpband)];

    % Make sure that cache data folder exists for this script
    [~, scriptName, ~] = fileparts(mfilename('fullpath'));
    scriptCacheFullpath = assertCacheDir(cfg.cacheDirFullpath, scriptName);
    fullCachefn = [scriptCacheFullpath matfnStr '.mat'];

    % Check if this run has been cached before
    isCached = checkCached(fullCachefn);
    
    % if current run of this function has a cached result already & not
    % overwriting
    if isCached && ~cfg.overwriteCache 
        load(fullCachefn);
        return % go ahead and skip all calculations and output the saved data
        
    end % if result is not cached, proceed with calculation below... 
 
end



%% Calculate running RMS on lfp data

% Get lfps
[lfp, timestamps, fs] = nwbLoadLfp(nwb, ...
    'channel', ch, ...
    'epoch', epochStr);


% Filter for band of interest
if isempty(lfpband) % default case, no extra band-filtering
    lfpFilt = lfp;

else
    [blfp, alfp] = butter(2, lfpband / (fs/2), 'bandpass');
    lfpFilt = bandfilt_edgecomp(blfp, alfp, double(lfp));

end


winSamp = floor(rmswin * fs);

% Get RMS over time
disp('BEGIN RMS SLIDING FILTER...')
tic
lfpRMS = smoothrms(lfpFilt, winSamp);
disp('DONE!!!')
toc
    
% Save results in cache if caching is turned on
if cfg.allowCache 
    save(fullCachefn, 'lfpRMS', 'timestamps', 'fs');
    
end


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
ch = namepars.ch;
LFP_BAND = namepars.LFP_BAND;
intPath = namepars.intPath;
% cpair = namepars.cpair;
WINSEC = namepars.WINSEC;

chStr = buildChStr(ch);

[~,nwblab,~] = fileparts(nwbfn);

formatSpecFn = 'RMS_%s_AvWin%ss_ch%s_lfpHz%s';
matfnStr = sprintf(formatSpecFn, ...
    nwblab, ...
    num2str(WINSEC), ...
    chStr, ...
    [num2str(LFP_BAND(1)) '-' num2str(LFP_BAND(2))]);

fullPathFn = [intPath matfnStr '.mat'];

end

function chStr = buildChStr(ch)
% converts numerical contents of ch to str, concats if more than one num

chStr = '';
for i = 1:length(ch)
    chStr = [chStr num2str(ch(i))];
    
end
    
end

function scriptCacheFullpath = assertCacheDir(cacheDirFullpath, scriptName)
% Get the name of the currently running script

% global PROJROOTPATH

% [~, scriptName] = fileparts(mfilename('fullpath'));

% First make sure that this script has a folder within the project folder's
% intermediate data section
scriptCacheFullpath = [cacheDirFullpath scriptName];
if ~exist(scriptCacheFullpath, 'dir')
    mkdir(scriptCacheFullpath)
    
end

scriptCacheFullpath = [scriptCacheFullpath '\'];

end

function xrms = smoothrms(x, winSamps)
% First get the square of every point, then smooth with a moving average
% square window, finally get the square root of every point. Note: x must
% be column vectors of data.

w = ones(winSamps,1);
w = w ./ sum(w);

xrms = x;
for i = 1:size(x, 2)
    xrms(:,i) = xrms(:,1).^2;
    xrms(:,1) = conv(xrms(:,1), w, 'same');
    xrms(:,1) = sqrt(xrms(:,1));

end

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
