% script for saving a downsampled version of all data within Stream Store
% of TDT block


%%


% remove TDTMatlabSDK matnwb version, and instead include my local version
% This is because matnwb current version has a bug when I try to read the
% harmaline tdt blocks, has to do with me having a TDT note file in the
% block

% This code assumes that NMRCNWBStandardization (and subdirectories) has 
% already been added to matlab path
% 
% addpath(genpath('C:\Users\bello043\Documents\GitHub\NMRCNWBStandardization'));
% 
% warning('off','all');
% rmpath(genpath('C:\Users\bello043\Documents\GitHub\NMRCNWBStandardization\toolbox\TDTMatlabSDK'));
% warning('on','all');
% addpath(genpath('C:\Users\bello043\Documents\GitHub\TDTMatlabSDK'));
% 
% 
% % Add project folder and all subdirectories within:
% addpath(genpath('C:\Users\bello043\Documents\GitHub\Thalamic-DBS-HarmalineTremor'));
% 
% % Remove .git portion of the file (some .git folders are huge/not needed
% % here), and suppress the MATLAB warning if the project folder is not
% % git-source-controlled:
% % warning('off', 'MATLAB:rmpath:DirNotFound');
% % rmpath([projRootPath, '\.git']);
% 
% % Add path to TDTmatlabSDK, which contains TDTbin2mat (loading TDT data)
% tdtmatlabsdkPath = 'C:\Users\bello043\Documents\GitHub\TDTMatlabSDK';
% addpath(genpath(tdtmatlabsdkPath));
% rmpath([tdtmatlabsdkPath, '\.git']);


% Add path to Git version-controlled code for Harmaline project
addpath(genpath('C:\Users\bello043\Documents\GitHub\Thalamic-DBS-HarmalineTremor'));

%% Load in table with metadata on all TDT blocks

pn = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Processing\LFPdownsamp\';
fn = 'TDTblockExpMetadata.xlsx';
T = readtable([pn fn]);



%% Get downsampled differential LFP

DATALOCPATH = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';

% save info
projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';
procDir = 'Data Processing\LFPdownsamp\';

storeStr = 'RAW8';


% for each TDT block, convert RAW8 to differential LFPs downsampled
nBlks = height(T);
for iBlk = 1:nBlks
    TDTTANK = [T.TDTtank{iBlk} '\']; 
    TDTBLOCK = T.TDTblock{iBlk}; 
    fulltdtpath = [DATALOCPATH TDTTANK TDTBLOCK];
    
    fullsavepath = [projRootPath procDir TDTBLOCK '_diffLFP'];
    
    
    % for each channel pair, load two raw channels, get Diff
    lfp = [];
    nChans = 8;
    for iCh = 1:(nChans-1)
        disp(['getting ' num2str(iCh) '-' num2str(iCh+1)])
        tic
        tdt = TDTbin2mat_chanDiff(fulltdtpath, [iCh, iCh+1], 'STORE', storeStr);
        rawdiff = tdt.streams.(storeStr).data';
        fs = tdt.streams.(storeStr).fs;

        % filt, then downsample
        fc = [0.5, 300];
        [b, a] = butter(2, fc / (fs/2), 'bandpass');
        filtdiff = filtfilt(b, a, double(rawdiff));

        factor = 48; % downsample by factor of 48 to bring fs to 1.0173e+03
        dwndiff = downsample(filtdiff, factor);
        dwnfs = fs / factor;

        % concat into final multichannel lfp matrix
        lfp = [lfp, dwndiff];

        toc

    end
    

    fs = dwnfs;
    save(fullsavepath, 'lfp', 'fs');

end


