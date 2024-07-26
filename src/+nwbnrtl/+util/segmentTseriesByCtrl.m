function [accSegFinal, tstSegFinal] = segmentTseriesByCtrl(acc, tst, ctrl, maxSamp, subSamp)
% Summary of this function goes here
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
    
    if ~exist('subSamp', 'var'), subSamp = maxSamp; end

    for iSeg = 1:nSegs
        if nSamps(iSeg) > maxSamp % Too big! subdivide the segment into smaller ones
            nSubs = floor(nSamps(iSeg)/subSamp) + 1;
            tempAcc = accSeg{iSeg,1};
            tempTst = tstSeg{iSeg,1};
            
            % mark subdivisions with indices
            xb = 1:(nSubs-1);
            bIdx = [1, xb*subSamp + 1];
            eIdx = [xb*subSamp, nSamps(iSeg)];
            
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



end

