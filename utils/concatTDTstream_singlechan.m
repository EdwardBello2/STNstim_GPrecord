% Concatenates TDT block stream data from several TDTblocks, outputs them
% in same struct format while ignoring any other data that the blocks may
% have, with output "tdtConcat" in the same struct format as the output
% of "TDTbin2mat.m". works with one channel at a time, of one specific
% stream field (e.g., tdt.streams.RAW8(1,:) ). Manipulates TDT block data in
% the struct format output by TDTbin2mat.m (see TDT matlab SDK for
% details). Data for concatenation must all have the same sampling
% frequency.
%
%
% SYNTAX:
% [tdtConcat, idxBlk, infoBlk] = concatTDTstream_singlechan(blockpath, tdtblocks, streamField, chan)
%
% 
% DESCRIPTION:
% [tdtConcat, idxBlk, infoBlk] = concatTDTstream_singlechan(blockpath, tdtblocks, streamField, chan)
% Specify the fullpath of where all TDT blocks of interest are currently
% located (must be in same folder, ideally the original TANK folder),
% Specify cell array of strings with filenames of blocks, and using char
% strings specify the field and channel of interest. Outputs the requested
% data common to all blocks in concatenated form. Discontinuities in data
% due to empty non-recorded time between blocks is not represented, so
% that sudden discontinuities will be present where concatenation occurs.
% Using the first block as time reference, "idxBlk" tracks 
% time-synchronized index position of each data point relative to first
% data point in first block. 
%
%
% INPUTS:
%   blockpath - fullpath to directory holding all TDT blocks, char array
%   tdtblocks - cell array of char arrays with names of all TDT blocks to
%               concat; they must be in order of occurrence. 
% streamField - name of stream, char array
%        chan - numerical index of channel of interest, single scalar value
% 
% 
% OUTPUTS: 
% tdtConcat - struct with one channel of stream data, so that data can be
%             accessed with dot-notation (e.g. tdtConcat.streams.data).
%    idxBlk - n x 2 array of first and last integer indices of all
%             re-synchronized data, for each of "n" blocks.  
%   infoBlk - n x 1 struct containing tdt.info fields from each of the "n"
%             blocks concatenated. 
% 
%
% EXAMPLES:
% 


% Author: Ed Bello
% Created: 3/26/2020
% Current ver: 0.2
%
% 
% ChangeLog:
% 0.2 -- changed it from filling inter-block space with NaNs to no fill,
% just simple concatenation, keeping track of the time-synchronized indices
% though
% 
%
% TO-DO

function [tdtConcat, idxBlk, idxBlkSync, infoBlk] = concatTDTstream_singlechan(BLOCKPATH, TDTBLOCKS, streamField, chan)
  
   
% storeField = 'RAW8'; 
assert(((numel(chan) == 1) && (chan > 0)), 'Only load one channel at a time!');

% subfunction to make sure blocks are appropriate to concatenate:
performCHECKS(BLOCKPATH, TDTBLOCKS, streamField);



nTDTBLKS = length(TDTBLOCKS);
idxBlk = zeros(nTDTBLKS, 2);
idxBlkSync = zeros(nTDTBLKS, 2);

% Initial values
dataConcat = [];
% concatSamps = 0;
% intNaNs = [];
disp(['BEGIN Concatenating ch' num2str(chan) ' of ' streamField '...'])
sampTot = 0;
for iBlk = 1:nTDTBLKS
%     fulltdtpath  = [DATALOCPATH, '\', TDTTANK, '\', TDTBLOCKS{iBlk}];
    fulltdtpath  = [BLOCKPATH, '\', TDTBLOCKS{iBlk}];

    
    % load tdt struct of current block
    tic
    tdt = TDTbin2mat(fulltdtpath, 'STORE', streamField, ...
    'CHANNEL', chan);
    toc
    dataBlk = tdt.streams.(streamField).data;
    
    % Get idxBlk for final concatenated data
    idxBlk(iBlk,1:2) = [1, length(dataBlk)] + sampTot;
    sampTot = idxBlk(iBlk,2);
    
    
    % Update the collection of tdt.info for each block
    infoBlk(iBlk) = tdt.info;
    
    
    % get initial start time and sasmpling freq of session-data
    if iBlk == 1
%         sess_tBegdt = datetime(tdt.info.utcStartTime);
        fs = tdt.streams.(streamField).fs; % samps/sec
        
    end  
    
    
    % obtain inter-block NaNs based on diff between this block and last
    if iBlk == 1
        intNaNs = [];
        
    else
        tBeg_currBlk = datetime(tdt.info.utcStartTime);
        intBlkSecs = seconds(tBeg_currBlk - tEnd_lastBlk);
        intBlkSamps = round(intBlkSecs * fs);
        intNaNs = single(NaN(1, intBlkSamps));
        
    end


    % add inter-block NaNs to dataSess
    dataConcat = [dataConcat, intNaNs];
    
    
    % update blkIdx mapping current block data to its position in datSess
    idxBlkSync(iBlk,1:2) = [1, length(dataBlk)]  + length(dataConcat);
    
    
    % add current block's data to dataSess
    dataConcat = [dataConcat, dataBlk];
    
    
    % get end-time of current block to serve in next iteration
    blk_secs = length(dataBlk) / fs; % seconds 
    tEnd_lastBlk = seconds(blk_secs) + datetime(tdt.info.utcStartTime);
       
end


% Remove NaN values from data, so that data has no inter-block filler
dataConcat(isnan(dataConcat)) = [];

% Create the output tdt struct as if it were a call to TDTbin2mat with all
% struct fields empty except "stream" where the new concatenated data from
% blocks is housed. 
tdtConcat = tdt; % just get the one from the most recent block
tdtConcat.epocs = [];
tdtConcat.info = [];
tdtConcat.streams.(streamField).data = dataConcat;


disp(['DONE Concatenating ch' num2str(chan) ' of ' streamField '!'])


end

%% SUB-FUNCTIONS

function TF = checkSTORE(BLOCKPATH, TDTBLOCKS, streamField)
% Check if all TDT blocks in question have "streamField" in common

nTDTBLKS = length(TDTBLOCKS);
hasStreamField = false(1, nTDTBLKS);

for iBlk = 1:nTDTBLKS
    fulltdtpath  = [BLOCKPATH, '\', TDTBLOCKS{iBlk}];
    tdt = TDTbin2mat(fulltdtpath, 'HEADERS', 1);
    hasStreamField(iBlk) = isfield(tdt.stores, streamField);
    
end

TF = all(hasStreamField);

end

function TF = checkSTOREFS(BLOCKPATH, TDTBLOCKS, streamField)
% Check if all TDT blocks in question have the same sampling frequency for
% streamField store

nTDTBLKS = length(TDTBLOCKS);
fsBlk = zeros(1, nTDTBLKS);

for iBlk = 1:nTDTBLKS
    fulltdtpath  = [BLOCKPATH, '\', TDTBLOCKS{iBlk}];
    tdt = TDTbin2mat(fulltdtpath, 'HEADERS', 1);
    fsBlk(iBlk) = tdt.stores.(streamField).fs;
    
end

TF = all(fsBlk == fsBlk(1));

end

function performCHECKS(BLOCKPATH, TDTBLOCKS, streamField)
% First make sure that all blocks do in fact have STORE field in common,
% under "streams"
streamPresent = checkSTORE(BLOCKPATH, TDTBLOCKS, streamField);
assert(streamPresent, ['Each TDT block must contain ', streamField]);
disp(' '); disp(['Confirmed: ', streamField, ' is present in all TDT blocks']); disp(' ');

% Make sure that stream STORE has the same sampling frequency
% across all blocks, otherwise concatenation is no bueno
sameFs = checkSTOREFS(BLOCKPATH, TDTBLOCKS, streamField);
assert(sameFs, [streamField, ' must have same sampling frequency across all blocks']);
disp(' '); disp(['Confirmed: ', streamField, ' sampling freq is consistent across blocks']); disp(' ');

end

