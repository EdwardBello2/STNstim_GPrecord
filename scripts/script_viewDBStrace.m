% script to look at subselection

clear; close all

iTDTblock = 'jdsktestv2-200124-143139.mat';
projRootPn = 'C:\DATAtemp\STNstim_GPrecord\';
acqPn = 'Data Acquisition\';
prcPn = 'Data Processing\';

% Load data tables with metadata for each block (BlockData.xlsx) and each 
% session composed of multiple blocks (SessionData.xlsx)
blkTab = readtable([projRootPn prcPn 'BlockData_Azula.xlsx']);
sessTab = readtable([projRootPn prcPn 'SessionData_Azula.xlsx']);
dataTab = join(blkTab, sessTab);

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
