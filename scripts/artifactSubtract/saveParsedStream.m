function saveParsedStream(matFullPath,streamStr,blockIntervals,saveDirFullPath)
% Convert data struct from TDTbin2mat into .mat with individual variables 
% for each channel in the store specified by "streamStr", and save to
% same folder path as matfn. 
% 
% Depends on "fieldnamesr.m" from Matlab Central file exchange to be on set path:
% Adam (2020). Get structure field names in recursive manner (https://www.mathworks.com/matlabcentral/fileexchange/33262-get-structure-field-names-in-recursive-manner), MATLAB Central File Exchange. Retrieved March 6, 2020.

% TO-DO:
% - Change the ^SUNx input of the save command at bottom to match whatever streamStr is 



% matfn = 'D:\PROJECTS\STNstim_GPrecord\Data Acquisition\012420mat\jdsktestv2-200124-160501.mat';
% streamStr = 'SUNx';
dataSubStructStr = [streamStr '.data'];
  fsSubStructStr = [streamStr '.fs'];
       ncharData = numel(dataSubStructStr);
         ncharFs = numel(fsSubStructStr);
       tdtstruct = load(matFullPath);

% Retrieve the data in the field specified in "streamStr" regardless of
% what the user decided to name the TDTbin2mat resultant struct (i.e. "data")
 NAMES = fieldnamesr(tdtstruct); % THIS FUNCTION MUST BE IN MATLAB PATH
nNames = length(NAMES);
isStructPathData = false(nNames, 1);
  isStructPathFs = false(nNames, 1);
for iName = 1:nNames
    iStr = NAMES{iName};
    try
        isStructPathData(iName) = strcmp(dataSubStructStr, iStr(end-ncharData+1:end));
        isStructPathFs(iName) = strcmp(fsSubStructStr, iStr(end-ncharFs+1:end));

    catch ME
        
    end           
    
end

FullStructStrData = NAMES{isStructPathData};
FullStructStrFs = NAMES{isStructPathFs};


% define streamdata and sampling frequency
eval(['streamData = tdtstruct.' FullStructStrData ';']);
eval(['streamfs = tdtstruct.' FullStructStrFs ';']);


% Separate streamdata into separate blocks based on user-defined time
% intervals, and save each block with individual variables for each channel:
[~, matfn, ~] = fileparts(matFullPath);

blkIdx = blockIntervals;
nBlks = size(blkIdx, 1)
for iBlk = 1:nBlks
    blkData = streamData(:,blkIdx(iBlk,1):blkIdx(iBlk,2));
    
    
    % Separate tdt multichannel streamdata into individual time series,
    % iteratively named:
    nChs = size(blkData, 1);
    for iCh = 1:nChs
        eval([streamStr '_ch' num2str(iCh) ' = blkData(' num2str(iCh) ',:);']);

    end

    saveFullPath = [saveDirFullPath matfn '_' streamStr '_blk' num2str(iBlk) '_chAll'];

    save(saveFullPath, '-regexp', ['^' streamStr]) 
    
    
end

disp('Done saving into SARGE-compliant format!');
disp(' ')
disp(['Original matfile: ' matfn ])
disp(['Extracted stream: ' streamStr])
disp(['num blocks: ' num2str(nBlks)])
disp(' ')




end