% script for processing IMU data stored within csv files into mat files
% in which all IMU data for a given session is contained (in one file
% rather than multiple), assuring that no negative time skips in time
% stamps occurred. 
%
% Final incorporation into NWB structure will require
% synchronization with the TDT epoc "RaspPi" which detects the moment at
% which the Raspberry Pi code started and began logging IMU data and
% timestamps locally  within itself. 

% Separate sections for dealing with all 8 harmaline experiment sessions

clear; 
dataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata';

projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';
savePath = [projRootPath 'Data Acquisition\Uva dose study\IMUdata\'];

% NOTE some of these origignal files are harder to deal with than others,
% due to various complications on recording days. Some have to be handled
% specially..



%% TremoLfpDBS-190927-100155, 2 mg/kg, expos 1
clear tdtBlocks tdtStructs IMUcell idxNegdt IMUsyncPulse filename

tdtTank = 'Uva-190927';
tdtBlock = 'TremoLfpDBS-190927-100155';
filename{1} = [tdtBlock, '_IMU.csv'];

% TDT blocks from that day that *could* have a RaspPi sync pulse within
tdtBlocks = {'TremoLfpDBS-190927-100155'};
tdt = TDTbin2mat([dataPath, '\', tdtTank, '\', tdtBlocks{1}], ...
    'TYPE', {'epocs'});
tdtStructs{1} = tdt;

% Extract from tdt files the sync pulse time(s) from the Raspberry Pi ,
% referenced according to session start time (first TDTblock)
IMUsyncPulse = getIMUsyncPulse(tdtStructs);


ifile = 1;
% load IMU timestamps and data from original csv file as a matlab table
fullLoadPath = [dataPath, '\', tdtTank, '\', filename{ifile}];
T = loadIMUorig(fullLoadPath);

% add sync rereference to align IMU timestamps with NWB file as table
% userdata
T.Properties.UserData.sync = IMUsyncPulse(ifile);

% Correct for any negative time skips
[IMUcell{ifile}, idxNegdt{ifile}] = correctNegTime(T);
disp(['# of negative time skips detected/corrected: ' num2str(numel(idxNegdt{ifile}))])


% save "clean" result as matlab table (.mat)
save([savePath tdtBlock, '_IMU'], 'IMUcell');
disp('Saved result as matlab table!')



%% TremoLfpDBS-191004-100637, 4 mg/kg, expos 1
clear tdtBlocks tdtStructs IMUcell idxNegdt IMUsyncPulse filename

tdtTank = 'Uva-191004';
tdtBlock = 'TremoLfpDBS-191004-100637';
filename{1} = [tdtBlock, '_IMU.csv'];

% TDT blocks from that day that *could* have a RaspPi sync pulse within
tdtBlocks = {'TremoLfpDBS-191004-100637'};
tdt = TDTbin2mat([dataPath, '\', tdtTank, '\', tdtBlocks{1}], ...
    'TYPE', {'epocs'});
tdtStructs{1} = tdt;

% Extract from tdt files the sync pulse time(s) from the Raspberry Pi ,
% referenced according to session start time (first TDTblock)
IMUsyncPulse = getIMUsyncPulse(tdtStructs);


ifile = 1;
% load IMU timestamps and data from original csv file as a matlab table
fullLoadPath = [dataPath, '\', tdtTank, '\', filename{ifile}];
T = loadIMUorig(fullLoadPath);

% add sync rereference to align IMU timestamps with NWB file as table
% userdata
T.Properties.UserData.sync = IMUsyncPulse(ifile);

% Correct for any negative time skips
[IMUcell{ifile}, idxNegdt{ifile}] = correctNegTime(T);
disp(['# of negative time skips detected/corrected: ' num2str(numel(idxNegdt{ifile}))])


% save "clean" result as matlab table (.mat)
save([savePath tdtBlock, '_IMU'], 'IMUcell');
disp('Saved result as matlab table!')



%% TremoLfpDBS-191011-104322, 6 mg/kg, expos 1
clear tdtBlocks tdtStructs IMUcell idxNegdt IMUsyncPulse filename

tdtTank = 'Uva-191011';
tdtBlock = 'TremoLfpDBS-191011-104322';
filename{1} = [tdtBlock, '_IMU.csv'];

% TDT blocks from that day that *could* have a RaspPi sync pulse within
tdtBlocks = {'TremoLfpDBS-191011-104322'};
tdt = TDTbin2mat([dataPath, '\', tdtTank, '\', tdtBlocks{1}], ...
    'TYPE', {'epocs'});
tdtStructs{1} = tdt;

% Extract from tdt files the sync pulse time(s) from the Raspberry Pi ,
% referenced according to session start time (first TDTblock)
IMUsyncPulse = getIMUsyncPulse(tdtStructs);


ifile = 1;
% load IMU timestamps and data from original csv file as a matlab table
fullLoadPath = [dataPath, '\', tdtTank, '\', filename{ifile}];
T = loadIMUorig(fullLoadPath);

% add sync rereference to align IMU timestamps with NWB file as table
% userdata
T.Properties.UserData.sync = IMUsyncPulse(ifile);

% Correct for any negative time skips
[IMUcell{ifile}, idxNegdt{ifile}] = correctNegTime(T);
disp(['# of negative time skips detected/corrected: ' num2str(numel(idxNegdt{ifile}))])


% save "clean" result as matlab table (.mat)
save([savePath tdtBlock, '_IMU'], 'IMUcell');
disp('Saved result as matlab table!')



%% TremoLfpDBS-191018-100615, 8 mg/kg, expos 1
clear tdtBlocks tdtStructs IMUcell idxNegdt IMUsyncPulse filename

tdtTank = 'Uva-191018';
tdtBlock = 'TremoLfpDBS-191018-100615';
filename{1} = [tdtBlock, '_IMU_1.csv'];
filename{2} = [tdtBlock, '_IMU_2.csv'];

% TDT blocks from that day that *could* have a RaspPi sync pulse within
tdtBlocks = {'TremoLfpDBS-191018-100615', ...
    'TremoLfpDBS-191018-123049', ...
    'TremoLfpDBS-191018-130103', ...
    'TremoLfpDBS-191018-133150', ...
    'TremoLfpDBS-191018-140250', ...
    'TremoLfpDBS-191018-142724'};
nBlocks = numel(tdtBlocks);
tdtStructs = cell(nBlocks, 1);
for iBlk = 1:nBlocks
    tdt = TDTbin2mat([dataPath, '\', tdtTank, '\', tdtBlocks{iBlk}], ...
    'TYPE', {'epocs'});
    tdtStructs{iBlk,1} = tdt;
    
end

% Extract from tdt files the sync pulse time(s) from the Raspberry Pi ,
% referenced according to session start time (first TDTblock)
IMUsyncPulse = getIMUsyncPulse(tdtStructs);


nFiles = numel(filename);
for ifile = 1:nFiles

    % load IMU timestamps and data from original csv file as a matlab table
    fullLoadPath = [dataPath, '\', tdtTank, '\', filename{ifile}];
    T = loadIMUorig(fullLoadPath);
    
    % add sync rereference to align IMU timestamps with NWB file as table
    % userdata
    T.Properties.UserData.sync = IMUsyncPulse(ifile);

    % Correct for any negative time skips
    [IMUcell{ifile}, idxNegdt{ifile}] = correctNegTime(T);
    disp(['# of negative time skips detected/corrected: ' num2str(numel(idxNegdt{ifile}))])

end


% save "clean" result as matlab table (.mat)
save([savePath tdtBlock, '_IMU'], 'IMUcell');
disp('Saved result as matlab table!')



%% TremoLfpDBS-191025-104651, 2 mg/kg, expos 2
clear tdtBlocks tdtStructs IMUcell idxNegdt IMUsyncPulse filename

tdtTank = 'Uva-191025';
tdtBlock = 'TremoLfpDBS-191025-104651';
filename{1} = [tdtBlock, '_IMU.csv'];

% TDT blocks from that day that *could* have a RaspPi sync pulse within
tdtBlocks = {'TremoLfpDBS-191025-104651'};
tdt = TDTbin2mat([dataPath, '\', tdtTank, '\', tdtBlocks{1}], ...
    'TYPE', {'epocs'});
tdtStructs{1} = tdt;

% Extract from tdt files the sync pulse time(s) from the Raspberry Pi ,
% referenced according to session start time (first TDTblock)
IMUsyncPulse = getIMUsyncPulse(tdtStructs);


ifile = 1;
% load IMU timestamps and data from original csv file as a matlab table
fullLoadPath = [dataPath, '\', tdtTank, '\', filename{ifile}];
T = loadIMUorig(fullLoadPath);

% add sync rereference to align IMU timestamps with NWB file as table
% userdata
T.Properties.UserData.sync = IMUsyncPulse(ifile);

% Correct for any negative time skips
[IMUcell{ifile}, idxNegdt{ifile}] = correctNegTime(T);
disp(['# of negative time skips detected/corrected: ' num2str(numel(idxNegdt{ifile}))])


% save "clean" result as matlab table (.mat)
save([savePath tdtBlock, '_IMU'], 'IMUcell');
disp('Saved result as matlab table!')



%% TremoLfpDBS-191101-101430, 4 mg/kg, expos 2
clear tdtBlocks tdtStructs IMUcell idxNegdt IMUsyncPulse filename

tdtTank = 'Uva-191101';
tdtBlock = 'TremoLfpDBS-191101-101430';
filename{1} = [tdtBlock, '_IMU.csv'];

% TDT blocks from that day that *could* have a RaspPi sync pulse within
tdtBlocks = {'TremoLfpDBS-191101-101430'};
tdt = TDTbin2mat([dataPath, '\', tdtTank, '\', tdtBlocks{1}], ...
    'TYPE', {'epocs'});
tdtStructs{1} = tdt;

% Extract from tdt files the sync pulse time(s) from the Raspberry Pi ,
% referenced according to session start time (first TDTblock)
IMUsyncPulse = getIMUsyncPulse(tdtStructs);


ifile = 1;
% load IMU timestamps and data from original csv file as a matlab table
fullLoadPath = [dataPath, '\', tdtTank, '\', filename{ifile}];
T = loadIMUorig(fullLoadPath);

% add sync rereference to align IMU timestamps with NWB file as table
% userdata
T.Properties.UserData.sync = IMUsyncPulse(ifile);

% Correct for any negative time skips
[IMUcell{ifile}, idxNegdt{ifile}] = correctNegTime(T);
disp(['# of negative time skips detected/corrected: ' num2str(numel(idxNegdt{ifile}))])


% save "clean" result as matlab table (.mat)
save([savePath tdtBlock, '_IMU'], 'IMUcell');
disp('Saved result as matlab table!')



%% TremoLfpDBS-191108-101829, 6 mg/kg, expos 2
clear tdtBlocks tdtStructs IMUcell idxNegdt IMUsyncPulse filename

tdtTank = 'Uva-191108';
tdtBlock = 'TremoLfpDBS-191108-101829';
filename{1} = [tdtBlock, '_IMU_1.csv'];
filename{2} = [tdtBlock, '_IMU_2.csv'];

% TDT blocks from that day that *could* have a RaspPi sync pulse within
tdtBlocks = {'TremoLfpDBS-191108-101829', ...
    'TremoLfpDBS-191108-124442', ...
    'TremoLfpDBS-191108-131346', ...
    'TremoLfpDBS-191108-134457', ...
    'TremoLfpDBS-191108-141530', ...
    'TremoLfpDBS-191108-144621'};
nBlocks = numel(tdtBlocks);
tdtStructs = cell(nBlocks, 1);
for iBlk = 1:nBlocks
    tdt = TDTbin2mat([dataPath, '\', tdtTank, '\', tdtBlocks{iBlk}], ...
    'TYPE', {'epocs'});
    tdtStructs{iBlk,1} = tdt;
    
end

% Extract from tdt files the sync pulse time(s) from the Raspberry Pi ,
% referenced according to session start time (first TDTblock)
IMUsyncPulse = getIMUsyncPulse(tdtStructs);


nFiles = numel(filename);
for ifile = 1:nFiles

    % load IMU timestamps and data from original csv file as a matlab table
    fullLoadPath = [dataPath, '\', tdtTank, '\', filename{ifile}];
    T = loadIMUorig(fullLoadPath);
    
    % add sync rereference to align IMU timestamps with NWB file as table
    % userdata
    T.Properties.UserData.sync = IMUsyncPulse(ifile);

    % Correct for any negative time skips
    [IMUcell{ifile}, idxNegdt{ifile}] = correctNegTime(T);
    disp(['# of negative time skips detected/corrected: ' num2str(numel(idxNegdt{ifile}))])

end


% save "clean" result as matlab table (.mat)
save([savePath tdtBlock, '_IMU'], 'IMUcell');
disp('Saved result as matlab table!')



%% TremoLfpDBS-191115-100127, 8 mg/kg, expos 2
clear tdtBlocks tdtStructs IMUcell idxNegdt IMUsyncPulse filename

tdtTank = 'Uva-191115';
tdtBlock = 'TremoLfpDBS-191115-100127';
filename{1} = [tdtBlock, '_IMU.csv'];

% TDT blocks from that day that *could* have a RaspPi sync pulse within
tdtBlocks = {'TremoLfpDBS-191115-100127'};
tdt = TDTbin2mat([dataPath, '\', tdtTank, '\', tdtBlocks{1}], ...
    'TYPE', {'epocs'});
tdtStructs{1} = tdt;

% Extract from tdt files the sync pulse time(s) from the Raspberry Pi ,
% referenced according to session start time (first TDTblock)
IMUsyncPulse = getIMUsyncPulse(tdtStructs);


ifile = 1;
% load IMU timestamps and data from original csv file as a matlab table
fullLoadPath = [dataPath, '\', tdtTank, '\', filename{ifile}];
T = loadIMUorig(fullLoadPath);

% add sync rereference to align IMU timestamps with NWB file as table
% userdata
T.Properties.UserData.sync = IMUsyncPulse(ifile);

% Correct for any negative time skips
[IMUcell{ifile}, idxNegdt{ifile}] = correctNegTime(T);
disp(['# of negative time skips detected/corrected: ' num2str(numel(idxNegdt{ifile}))])


% save "clean" result as matlab table (.mat)
save([savePath tdtBlock, '_IMU'], 'IMUcell');
disp('Saved result as matlab table!')



%% SUB-FUNCTIONS

function IMUsyncPulse = getIMUsyncPulse(tdtStructs)
% takes as input a cell array of tdt data structures extracted from all the
% day's TDT blocks using TDTbin2mat, in chronological order of recording

nBlocks = numel(tdtStructs);
iBlock = 1;

for iBlock = 1:nBlocks
    tdt = tdtStructs{iBlock,1};
    utcStartTime{iBlock} = tdt.info.utcStartTime;

end

% Determine how to re-reference times according to first TDT block

    
% CREATE timestamps vector that is based on global syncronized reference
% time. If data has missing chunks of time, the time-skips will be
% reflected by discontinuities in the timestamps, according to NWB standard
% format guidelines. 

refDT = datetime(utcStartTime(1));
IMUsyncPulse = [];
for iBlock = 1:nBlocks
    tdt = tdtStructs{iBlock,1};
    
    if isfield( tdt.epocs, 'RsP_')
        % Get non-synchronized RaspPi sync time for each TDT block
        sync = tdt.epocs.RsP_.onset;

        % Get synchronization reference time (seconds) for each TDT block
        refSec = seconds(datetime(utcStartTime(iBlock)) - refDT);

        % Re-reference4 all timestamps accordingly
        syncRef = sync + refSec;

        IMUsyncPulse = [IMUsyncPulse; syncRef];
    
    end

end




end

function [IMU, idxNegdt] = correctNegTime(T)
% Timestamps should have monotonically increasing time, however at least
% one file had an issue where time restarted ~0 all of a sudden, thus there
% was a negative skip in time. This function compensates for such jumps and
% returns time to monotonically increasing. Output is table identical to 
% input if no negative timeskips occur. If there are any timeskips, then
% output is mostly identical to original, except: 1) timestamps now all 
% rereferenced, and 2) each time one skip occured, we remove one timestamp
% along with all corresponding data points. 

IMU = T;
tst = IMU.TimeFromStart_s_; 


% testdata = 1:numel(tst);
% first detect any negative time skips
dt = diff(tst);
% isNeg = dt < 0;
idxNegdt = find(dt < 0);

% correct timestamps to match back up and cancel any negative skips
looptst = tst;

nIdx = numel(idxNegdt);
for i = 1:nIdx
    currentTst = looptst(idxNegdt(i));
    looptst(idxNegdt(i)+1) = 0;
    looptst((idxNegdt(i)+1):end) = looptst((idxNegdt(i)+1):end) + currentTst;
    
end


% remove duplicate timestamps and corresponding data points
% looptst(idxNegdt) = [];
IMU.TimeFromStart_s_(:) = looptst;
IMU(idxNegdt,:) = [];



end

function [T] = loadIMUorig(fullLoadPath)
% Detect current import options and remove the first column from import
% (col 1 is a datetime object that's not useful here)
opts = detectImportOptions(fullLoadPath);
opts.SelectedVariableNames(1) = [];


T = readtable(fullLoadPath, opts);
disp('done reading in table')

end

