function [ScleanCell, tstSegHarRef, fSegs, resultsInfo] = func_gatherACCproc_allPSDTimeSegs_har(nwb, cfg)
% function that gathers artifact-clean PSDs for each recorded file based on
% spectrogram time windows. Large files are broken up into smaller segments
% (cfg fields). Artifact-filled spectrogram windows are identified (by detecting massive
% spectral leakage) and are removed from the collection within a segment of data.
% If specified a segment may be discarded if too short or if too few useful
% spectrogram windows (i.e. individual PSDs in time) are left (cfg fields).
%

% 
% myStatFun = @(x)mean(x);


procModule = nwb.processing.get('ACCproc').deref();
accTseriesObj = procModule.nwbdatainterface.get('ACCproc');

% accTseriesObj = nwb.processing.get('ACCproc').deref();


% load just the accelerometry data
dims = accTseriesObj.data.dims;
acc = accTseriesObj.data.load([1, 1], [dims(1), 3]);
tst = accTseriesObj.timestamps.load;
ctrl = accTseriesObj.control.load;


fs = accTseriesObj.starting_time_rate;
% maxSegTime = 300; % seconds
tooBigSamp = floor(fs * cfg.tooBigTime);
subSamp = floor(fs * cfg.subSegTime);
[accSeg, tstSeg] = nwbnrtl.util.segmentTseriesByCtrl(acc, tst, ctrl, ...
    tooBigSamp, subSamp);
resultsInfo.numDataSegments.original = length(accSeg);

% Extra step to remove any segments that are too short in time
[accSeg, tstSeg, isTooShort] = remSegsTooShort(accSeg, tstSeg, cfg.tooSmallTime);
resultsInfo.numDataSegments.tooShort = sum(isTooShort);

nSegs = length(accSeg);


% Obtain summed triaxial spectrogram over time
for iSeg = 1:nSegs
    [SaccCell{iSeg,1}, fSegs{iSeg,1}, tSpec] = spectrogramacc3(accSeg{iSeg,1}, fs, cfg);
       
end


% Clean spectrograms of any movement artifact spectral leakage
% make it based on whole-recording rather than local segments
% for iSeg = 1:nSegs
[ScleanCell, resultsInfo.movArtifacts.RemovedBySegment ...
    resultsInfo.movArtifacts.artIdxSeg] = removeMovArtInAll(SaccCell, cfg);
    
% end

% f = fSegs{1};
% % Calculate tremor power for specified tremor band
% isTremor = (f >= cfg.tremorBand(1)) & (f < cfg.tremorBand(2));
% for iSeg = 1:nSegs
%     Stemp = ScleanCell{iSeg,1};
%     tremPow{iSeg,1} = sum(Stemp(isTremor,:));
%     
% end
% ScleanCell = ScleanCell;

% Adjust timestamps to reflect harmaline injection time (t = 0)
dstr = datestr(nwb.timestamps_reference_time);
startDT = datetime(dstr(end-7:end));

harStr = nwb.general_notes(end-7:end);
harDT = datetime(nwb.general_notes(end-7:end));
harRefSec = seconds(harDT - startDT);

for iSeg = 1:nSegs
    tstSegHarRef{iSeg,1} = tstSeg{iSeg} - harRefSec;
    
end

% Remove any segments of data that have too few fband power data points
[ScleanCell, tstSegHarRef, isTooSparse] = remSegsSparseData(ScleanCell, tstSegHarRef, cfg.tooSparseData);
resultsInfo.numDataSegments.tooSparse = sum(isTooSparse);

resultsInfo.numDataSegments.final = length(ScleanCell);


% provide for tracking data segments of original data, before removal of
% movement artifact windows from both acc data and time segments
resultsInfo.origDataSeg = SaccCell;
resultsInfo.origTstSeg = tstSeg;

end



%% SUB-FUNCTIONS

function [accSeg, tstSeg, isRemove] = remSegsSparseData(accSeg, tstSeg, tooSparseData)
% Mark the segments that have too little data points left after removing
% movement artifact data. then remove them from the cell array

nSegs = length(tstSeg);
isRemove = false(nSegs,1);
for iSeg = 1:nSegs
    acc = accSeg{iSeg};
    if size(acc, 2) < tooSparseData
        isRemove(iSeg) = true;
        
    end
    
end

accSeg(isRemove) = [];
tstSeg(isRemove) = [];


end

function [accSeg, tstSeg, isRemove] = remSegsTooShort(accSeg, tstSeg, tooSmallTime)
% Mark the segments that are too short in time width, then remove them from
% the cell array

nSegs = length(tstSeg);
isRemove = false(nSegs,1);
for iSeg = 1:nSegs
    tst = tstSeg{iSeg};
    if (tst(end) - tst(1)) <= tooSmallTime
        isRemove(iSeg) = true;
        
    end
    
end

accSeg(isRemove) = [];
tstSeg(isRemove) = [];



end

function [accSegFinal, tstSegFinal] = segmentacc(acc, tst, ctrl, maxSamp)
% separate acc data and timestamps into time segments based on recordings.
% Large recordings are further subdivided into smaller time segements according
% to max segment time. maxTime (seconds). ctrl has zero at the start of
% every recording. 

% for acc and tst composed of N separate recordings, ctrl tracks which data
% belongs to which N recording

% detect the number of segments
nSegs = max(ctrl);
for iSeg = 1:nSegs
    accSeg{iSeg,1} = acc((ctrl == iSeg),1:3);
    tstSeg{iSeg,1} = tst((ctrl == iSeg));
    nSamps(iSeg) = sum(ctrl == iSeg);
    
end

% If any segments are too long according to maxSamp 
accSegFinal = {};
tstSegFinal = {};
if exist('maxSamp', 'var')
    for iSeg = 1:nSegs
        if nSamps(iSeg) > maxSamp % Too big! subdivide the segment into smaller ones
            nSubs = floor(nSamps(iSeg)/maxSamp) + 1;
            tempAcc = accSeg{iSeg,1};
            tempTst = tstSeg{iSeg,1};
            
            % mark subdivisions with indices
            xb = 1:(nSubs-1);
            bIdx = [1, xb*maxSamp + 1];
            eIdx = [xb*maxSamp, nSamps(iSeg)];
            
            % loop-fill the new subsections to be added
            for iSub = 1:nSubs
                accSegFinal = [accSegFinal; {tempAcc(bIdx(iSub):eIdx(iSub),1:3)}];
                tstSegFinal = [tstSegFinal; {tempTst(bIdx(iSub):eIdx(iSub))}];
                
            end
            
%             xb = [1:nRowsAdd]';
% 
%             % begIdxCorr
%             begIdxSmall = begIdx(i) + xb*maxSamp + 1;
%             begIdxCorr = [begIdxCorr; begIdx(i); begIdxSmall];
% 
%             % endIdxCorr
%             endIdxSmall = begIdx(i) + xb*maxSamp;
%             endIdxCorr = [endIdxCorr; endIdxSmall; endIdx(i)];

        else
            accSegFinal = [accSegFinal; accSeg(iSeg,1)];
            tstSegFinal = [tstSegFinal; tstSeg(iSeg,1)];
        
        end
    
        
        
    end
    
else % simple case where maxSamp is not entered into the function, no subdivision done
    accSegFinal = accSeg;
    tstSegFinal = tstSeg;
    
end




% 
% 
% begIdx = find(ctrl == 0);
% endIdx = [(begIdx(2:end)-1); length(ctrl)];
% % nSamps = diff([begIdx endIdx], [], 2);
% 
% % enforce maximum samples for each segment
% % isToobig = nSamps > maxSamp;
% nSamps = diff([begIdx endIdx], [], 2);
% begIdxCorr = [];
% endIdxCorr = [];
% for i = 1:length(begIdx)
% 
%     if diff([begIdx(i) endIdx(i)], [], 2) > maxSamp % Too big! break up the segment into smaller ones
%         nRowsAdd = floor(nSamps(i)/maxSamp); 
%         xb = [1:nRowsAdd]';
%         
%         % begIdxCorr
%         begIdxSmall = begIdx(i) + xb*maxSamp + 1;
%         begIdxCorr = [begIdxCorr; begIdx(i); begIdxSmall];
%         
%         % endIdxCorr
%         endIdxSmall = begIdx(i) + xb*maxSamp;
%         endIdxCorr = [endIdxCorr; endIdxSmall; endIdx(i)];
%         
%     else
%         begIdxCorr = [begIdxCorr; begIdx(i)];
%         endIdxCorr = [endIdxCorr; endIdx(i)];
%         
%     end
%     
%     
%     
% end


% Gather data segments
% nSegs = length(begIdxCorr);
% for iSeg = 1:nSegs
%     accSeg{iSeg,1} = acc(begIdxCorr(iSeg):endIdxCorr(iSeg),:);
%     tstSeg{iSeg,1} = tst(begIdxCorr(iSeg):endIdxCorr(iSeg),:);
%     
% end


end

function [accProc, tstProc] = fillACC(acc, tst, ctrl, fs)
% pre-process acc data so that missing values are interpolated and
% inter-file spaces are NaNs spaced out according to common
% sampling-frequency of all files in session
% ctrl goes from 0 to N for a single file of N samples (some may be
% missing). For an NWB file where multiple files were concatenated, the 0
% to N ctrl count will reset and start again from 0 at the startint point
% of the next file. 


Ts = 1/fs;                           % Sampling Time
% tst = [0:10  14:20  25:30]*0.00274;       % Create Time Series With Gaps
% acc = rand(size(tst));                      % Create Matching Data
tst = tst';
tn = round(tst/Ts);                       % Create Indices With Gaps
dt = diff([0 tn]);                      % Create Vector Of Differences
tg = find(dt > 1);                      % Find Indices Of Gaps
gaps = dt(tg)-1;                        % Find Lengths Of Gaps
ti = linspace(min(tst),max(tst),max(tn));   % Create Continuous Time Vector
fi = interp1(tst,acc,ti);                   % Create Matching Data Vector
for k1 = 1:length(tg)                   % Loop Filling Gaps In ‘f’ With NaN
    q = [tg(k1):tg(k1)+gaps(k1)-1];
    fi(tg(k1):tg(k1)+gaps(k1)-2) = nan(1,gaps(k1)-1);
end


% First subdivide 




nCtrls = length(ctrl);
for ictrl = 1:nCtrls
    
    
end

% then interpolate any short missing values within recordings


end

function [Sacc, f, tSpec] = spectrogramacc3(acc, fs, cfg)

% Detrend data, linear best fit is subtracted from each channel
acc = detrend(acc);


% Highpass filter out slow components

[b,a] = butter(3, cfg.fcDetrend/(fs/2), 'high');
accFilt = filtfilt(b, a, acc);



%% Get combined triaxial spectrogram

window = floor(cfg.TWIN * fs);

% Sum PSDs from all 3 accel channels
clear Z S
for ch = 1:3
    [Z(:,:,ch), f, tSpec] = spectrogram(accFilt(:,ch), window, cfg.NOVERLAP, [], fs, 'yaxis');
    % tRef = linspace(tSpec(1), tSpec(end), numel(tSpec));
    S(:,:,ch) = abs(Z(:,:,ch));
    S(:,:,ch) = S(:,:,ch).^2;

end
% Combine spectral power from all three axes (summate)
Sacc = sum(S, 3);

end

function [ScleanCell, nRemovedBySegment, isMovArtCell] = removeMovArtInAll(SaccCell, cfg)
% First concatenate all segments spectrograms together and keep track of
% their place within segment cells
nSegs = length(SaccCell);
SaccAll = [];
origIdx = [];
for iSeg = 1:nSegs
    iSacc = SaccCell{iSeg,1};
    origIdx = [origIdx, 1:size(iSacc, 2)]; % keeps tabls of time windows in spectrograms
    SaccAll = [SaccAll, SaccCell{iSeg,1}];   
    
end


% Identify movement artifact windows in spectrogram

[~, isMovArt] = util.remoutliers(sum(SaccAll, 1), ...
    'bound', 'upper', ...
    'MADthresh', cfg.MADthresh);
temp = SaccAll;
temp(:,isMovArt) = NaN;
SaccMarkOutliers = SaccAll;
SaccMarkOutliers(:,isMovArt) = NaN;


% Show Spectrogram with detected movement artifacts marked as 100's values
if cfg.displayWorkflowFigs
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
    title('Mov Artifacts marked');
    ylabel('Frequency (Hz)');
    f4.Position = [2555 57 560 421];

end

% Show Spectrogram with movement artifacts removed totally
if cfg.displayWorkflowFigs
    f4 = figure;
    ax4 = axes;
    ax4.Parent = f4;
%     tset = t;
%     tset(isMovArt) = [];
    surf(1:size(temp, 2), f, 10*log10(temp), ...
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
    title('Mov Artifacts removed')
    ylabel('Frequency (Hz)');
    f4.Position = [2554 566 560 420];

end


% Get average PSD free of movement artifact windows
ScleanAll = temp; % movement artifact windows removed
% PSD = mean(ScleanAll, 2);
if cfg.displayWorkflowFigs
    f1 = figure; ax = axes;
    plot(f, (PSD));
    grid on; 
    ax.XLim = [0, 35];
    %     ax.YLim = [0, popMaxPSD];
    xlabel('Frequency (Hz)');
%     title(['Minutes since harmaline inject: ' num2str(minutes(sessTabRow.harRefTime(1)))], 'Interpreter', 'none')

end


% Output nan-marked spectrograms back in their original cells
begIdx = find(origIdx == 1);
endIdx = [(begIdx(2:end)-1), size(ScleanAll, 2)];
for iSeg = 1:nSegs
    ScleanCell{iSeg,1} = ScleanAll(:,(begIdx(iSeg):endIdx(iSeg)));
    isMovArtCell{iSeg,1} = isMovArt(begIdx(iSeg):endIdx(iSeg));
    
end

% Remove nan-marked columns in each spectrogram (pertaining to time-window
% labeled as movemend artifact) for each segment
for iSeg = 1:nSegs
    iSaccnan = ScleanCell{iSeg,1};
    nanIdx = isnan(sum(iSaccnan, 1));
    iSaccClean = iSaccnan;
    iSaccClean(:,nanIdx) = [];
    
    ScleanCell{iSeg,1} = iSaccClean;
    nRemovedBySegment(iSeg,1) = sum(nanIdx);
    
end



end