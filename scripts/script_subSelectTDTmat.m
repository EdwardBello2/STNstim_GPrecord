% script to look at subselection

clear; clc; close all

iTDTblock = 'jdsktestv2-200203-123234.mat';
projRootPn = 'C:\DATAtemp\STNstim_GPrecord\';
acqPn = 'Data Acquisition\';
prcPn = 'Data Processing\';

% Load data tables with metadata for each block (BlockData.xlsx) and each 
% session composed of multiple blocks (SessionData.xlsx)
blkTab = readtable([projRootPn prcPn 'BlockData_Azula.xlsx']);
blkTab2 = readtable([projRootPn prcPn 'BlockData_Azula_250uA.xlsx']);
sessTab = readtable([projRootPn prcPn 'SessionData_Azula.xlsx']);
dataTab = join(blkTab, blkTab2, 'Keys', 'TDTblock');
dataTab = join(dataTab, sessTab);

iBlkTab = dataTab(strcmp(dataTab.TDTblock, iTDTblock),:);

spkChans = str2num(iBlkTab.spkChans{:}); % Specifies which channels to process
DBSelec = iBlkTab.DBSelec{:};


% load data and plot DBS detection channel
TDTblock = iBlkTab.TDTblock{:};
TDTtank = [iBlkTab.TDTtank{:} '\'];
load([projRootPn acqPn TDTtank TDTblock], 'data');

ch = iBlkTab.stimDetectChan;
store = iBlkTab.stimDetectStore{:};
dbsTrace = data.streams.(store).data(ch,:);


fs = data.streams.(store).fs;
tStamp = (1/fs) * (0:(length(dbsTrace) - 1));



figure; plot(tStamp, dbsTrace)

%% Create a sub-selected version of the data in the original mat file, save

% Specify the time limits to subselect data
tLim = [250, 330];

% iterate thru each store to subselect
names = fieldnames(data.streams);
nStores = length(names);
for iStore = 1:nStores
    data2.streams.(names{iStore}) = data.streams.(names{iStore}); % init new tdt file portion as copy
    
    idata = data.streams.(names{iStore}).data;
    fs = data.streams.(names{iStore}).fs;
    tStamp = (1/fs) * (0:(length(idata) - 1));
    
    idxSub = (tStamp >= tLim(1)) & (tStamp < tLim(2));
    
    idataSub = idata(:,idxSub);
    data2.streams.(names{iStore}).data = idataSub;
    
    
end

data2.info = data.info;
data2.time_ranges = data.time_ranges;

data = data2;
save([projRootPn acqPn TDTtank TDTblock(1:end-4) '_250uA'], 'data');
disp(['COMPLETE: ' TDTblock(1:end-4) '_250uA']);






