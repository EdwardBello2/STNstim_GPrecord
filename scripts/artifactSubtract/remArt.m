function [dataClean] = remArt(data, fs, stimSamps, params)
% for now I have Global Av and Dynamic Av specified as possible artifact
% removal strategies...


% nTempl = 30;

% blankSamps = 20; % samples to the right of stimSamp time to blank
blankSamps = params.blankSamps;

medSeg = median(diff(stimSamps));

% get segments around stim times
[artSegs,segIdx] = segments(data, stimSamps, medSeg);


%--------------------------------------------------------------------------
% SUBTRACT ARTIFACTS:

% Global Average:
% subtrSegs = artSub_globalAv(artSegs);

% Dynamic Average:
kMean = 30;
[subtrSegs, meanSeg] = artSub_DynamicAv(artSegs, kMean);


%--------------------------------------------------------------------------


% Blank specified region around pulse
blankSegs = subtrSegs;
blankSegs(:,1:blankSamps) = 0;
% figure; plot(blankSegs')



%% Segment data into arbitrary-length segments

[segData, segDataIdx] = segmentCells(data, stimSamps);
    


%% Insert artifact-removed segments into Segmented data

nSegs = size(artSegs, 1); % same as for segData
for iSeg = 1:nSegs
    % Replace data in dataSeg with as much of template-sub data as we can
    maxSegLength = min(length(blankSegs(iSeg,:)), length(segData{iSeg}));
    dataRec = segData{iSeg};
    dataRec(1:maxSegLength) = blankSegs(iSeg,1:maxSegLength);
    
    
    % if there are any remaining points in dataRec not covered by the
    % template-subtraction, shift points up or down to line up with last
    % point of altered portion:
    if length(dataRec) > maxSegLength
        diffData = dataRec(maxSegLength) - dataRec(maxSegLength+1);
        dataRec((maxSegLength+1):end) = dataRec((maxSegLength+1):end) +  diffData;
        
    end
    
    
    segData{iSeg} = dataRec;
%     dataIdx = segIdx(iSeg,1) - 1;
%     diffData = dataReconstruct(dataIdx) - blankSegs(iSeg,1);
%     iSegAdj = blankSegs(iSeg,:) + diffData;
%     dataReconstruct(segIdx(iSeg,:)) = iSegAdj;
%     
%     
%     % at end of this segment, adjust all values in signal to continue 
%     dataIdx = segIdx(iSeg,end) + 1;
%     diffData = dataReconstruct(dataIdx-1) - dataReconstruct(dataIdx);
%     dataReconstruct(dataIdx:end) = dataReconstruct(dataIdx:end) + diffData;

end




%% re-assemble segmented signal continuously 

dataReconstruct = data;
iSeg = 1;

% Make sure first segment's idx=1 data point equals the immediately
% preceding sample in original data
diffData = dataReconstruct(segDataIdx{iSeg}(1) - 1) - segData{iSeg}(1);

% replace data with processed segment, shifted by diffData
dataReconstruct(segDataIdx{iSeg}) = segData{iSeg} + diffData;


nSegs = size(artSegs, 1);
for iSeg = 2:nSegs
    % for each segment, compare first sample with preceding data sample, adjust
    % values in segment to be level
    diffData = dataReconstruct(segDataIdx{iSeg}(1) - 1) - segData{iSeg}(1);
    
    % replace data with processed segment, shifted by diffData
    dataReconstruct(segDataIdx{iSeg}) = segData{iSeg} + diffData;
    

end


% Append remaining data to the end of final processed segment, adjusting
% for height difference as before
diffData = dataReconstruct(segDataIdx{end}(end)) - ...
    dataReconstruct(segDataIdx{end}(end)+1);

dataReconstruct((segDataIdx{end}(end)+1):end) = ...
    dataReconstruct((segDataIdx{end}(end)+1):end) + diffData;



%%  final filter

fc = params.fc;
if numel(fc) == 2
    [b,a] = butter(2, fc / (fs/2), 'bandpass');
    
elseif numel(fc) == 1
    [b,a] = butter(2, fc / (fs/2), 'high');
    
else
    error('fc needs one or two values')
    
end

dataClean = filtfilt(b, a, dataReconstruct);
% plot(dataClean);

[clnSegs,segIdx] = segments(dataClean, stimSamps, medSeg);

figure; plot(clnSegs')


end % END of FUNCTION

%% SUB-FUNCTIONS


function [Segs,SegIndices] = segments(sig,refs,segLength)
% Organize a single-channel time-series into segments based off of
% periodic reference events. Note: window idicies must not exceed the
% indices of the original signal. 
%
% [Segs,SegIndices] = segments(sig,refs,segLength)
%       Output is an m-by-n matrix of segment windows according to "m"
%       number of reference events, as well as the indices of these
%       segments from the original time-series vector. Window-boundaries:
%       ref-event + "segLength" to the right.
% [Segs,SegIndices] = segments(sig,refs,[loffset, roffset])
%       Same as above, except that segment windows boundaries are defined
%       "loffset" samples to the left of ref-event and "roffset" samples to
%       the right of ref-event.

if ~isvector(sig), error('input "sig" must be a single vector value!'); end
if ~isvector(refs), error('input "refs" must be a single vector vaclue!'); end

%check segLength input
if length(segLength)>1
    if length(segLength)>2
        error('input exceeds allowed dimensions');
    end
    loffset = segLength(1);
    roffset = segLength(2);
elseif isscalar(segLength)
    loffset = 0;
    roffset = segLength-1;
else
    error('improper input for segLength')
end

refCount = length(refs);
Segs = zeros(refCount,loffset+1+roffset);
SegIndices = zeros(refCount,loffset+1+roffset);

for i = 1:refCount
    try
    SegIndices(i,:) = (refs(i)-loffset):(refs(i)+roffset);
    Segs(i,:) = sig(SegIndices(i,:));
    
    catch
        error(['i == ' num2str(i)])
        
    end
end  
% if nargout == 1
%     varargout{1} = SegIndices;
% end


end

function [segCells, segIdx] = segmentCells(sig, refs)
% Rather than have a matrix of equal-length segments like "segments.m", the
% idea here is to have a cell array of segment-ish data of arbitrary length
% from segment to segment, based only on the sample indices in "refs". 

if ~isvector(sig), error('input "sig" must be a single vector value!'); end
if ~isvector(refs), error('input "refs" must be a single vector vaclue!'); end


dSamp = median(diff(refs));

refCount = length(refs);
segCells = cell(refCount, 1);
segIdx = cell(refCount, 1);

for i = 1:(refCount-1)
    segIdx{i} = refs(i):refs(i+1);
    segCells{i} = sig(refs(i):refs(i+1));
    
end  

% Specify extent of segment after last ref-event:
if (refs(end) + dSamp) <= length(sig)
    segIdx{end} = refs(end):(refs(end) + dSamp);
    segCells{end} = sig(refs(end):(refs(end) + dSamp));
    
else
    segIdx{end} = refs(end):length(sig);
    segCells{end} = sig(refs(end):length(sig));   
    
end
% if nargout == 1
%     varargout{1} = SegIndices;
% end





end

function subtrSegs = artSub_globalAv(artSegs)

% Create Average template
artTempl = mean(artSegs,1);

% subtract from each segment
subtrSegs = artSegs - artTempl;

end

function [subtrSegs, meanSeg] = artSub_DynamicAv(artSegs, kMean)
% based closely on SARGE implementation of Dynamic average (check out their
% article for details, or take a peek in their code...

nStims = size(artSegs, 1);
meanSeg = artSegs;
subtrSegs = artSegs;

for iStim = 1:nStims
    pre = min(iStim - 1,kMean);
    post = min(nStims - iStim,kMean);
    meanSeg(iStim,:) = mean(artSegs(iStim-pre:iStim+post,:));
    subtrSegs(iStim,:) = artSegs(iStim,:) - meanSeg(iStim,:);

end

end