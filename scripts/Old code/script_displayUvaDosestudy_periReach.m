% script to quantify tremor-presence during TS task reaches


% Load nwb file
% Loading and saving pathways
clear nwb
nwbDataPath = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\NWBdata\';
projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';

% Uva dose study data of choice:
nwbFullPath = nwbDataPath;
% SESSION_ID = 'TremoLfpDBS-190927-100155'; % 1st 2 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191025-104651'; % 2nd 2 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191004-100637'; % 1st 4 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191101-101430'; % 2nd 4 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191011-104322'; % 1st 6 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191108-101829'; % 2nd 6 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191018-100615'; % 1st 8 mg/kg day
% SESSION_ID = 'TremoLfpDBS-191115-100127'; % 2nd 8 mg/kg day

% Baselien day
SESSION_ID = 'TremoLfpDBS-190923-110101';
% s = 8;


% Remove Chronux from Matlab's searchpath to prevent a problem with
% jacknife function used by bootci
rmpath(genpath('C:\Users\bello043\Documents\GitHub\Thalamic-DBS-HarmalineTremor\toolboxes\chronux_2_10'))


% SESSIONID_labels = {...
%     'TremoLfpDBS-190927-100155', ...
%     'TremoLfpDBS-191025-104651', ...
%     'TremoLfpDBS-191004-100637', ...
%     'TremoLfpDBS-191101-101430', ...
%     'TremoLfpDBS-191011-104322', ...
%     'TremoLfpDBS-191108-101829', ...
%     'TremoLfpDBS-191018-100615', ...
%     'TremoLfpDBS-191115-100127'};
% 
% SESSIONID_details = {...
%     '1st 2 mg/kg', ...
%     '2nd 2 mg/kg', ...
%     '1st 4 mg/kg', ...
%     '2nd 4 mg/kg', ...
%     '1st 6 mg/kg', ...
%     '2nd 6 mg/kg', ...
%     '1st 8 mg/kg', ...
%     '2nd 8 mg/kg'};


% Load nwb master file for this session, get acc TimeSeries object
disp('Reading NWB file...')
nwbM = nwbRead([nwbFullPath SESSION_ID '_MASTER.nwb']);
disp('DONE!')



% gather COMPLETED reaches timing data
behavEvents = nwbM.acquisition.get('PeriReachTimes').deref();

reachEndTseries = behavEvents.timeseries.get('reachEnd');
reachEnd = reachEndTseries.timestamps.load;
reachOnsetTseries = behavEvents.timeseries.get('reachOnset');
reachOnset = reachOnsetTseries.timestamps.load;
targetOnsetTseries = behavEvents.timeseries.get('targetOnset');
targetOnset = targetOnsetTseries.timestamps.load;


% Display peri-reach accel data
procModule = nwbM.processing.get('ACCproc').deref();
accTseriesObj = procModule.nwbdatainterface.get('ACCproc');

% load just the accelerometry data
dims = accTseriesObj.data.dims;
acc = accTseriesObj.data.load([1, 1], [dims(1), 3]);
tst = accTseriesObj.timestamps.load;
fs = accTseriesObj.starting_time_rate;
% ctrl = accTseriesObj.control.load;


% bandpass the accel data 
fc = [8 12];
[b,a] = butter(4, fc / (fs/2), 'bandpass');
accFilt = filtfilt(b, a, acc);

% reachWin = [-2 2];
% nReaches = length(reachEnd);
% for iReach = 1:nReaches
%     isPeriReach = (tst > (reachEnd(iReach) + reachWin(1))) & ...
%         (tst <= (reachEnd(iReach) + reachWin(2)));
%     figure; 
%     accReach = accFilt(isPeriReach,:);
%     t = (1/fs) * (0:(length(accReach)-1));
%     plot(t, accReach);
%     
% end


%% calc peri-reach tremor power

% for each reach
% reachWin = [-2 2];
accMag = abs(hilbert(accFilt));
% accMag = accFilt;


nReaches = length(reachEnd);
reachPow = zeros(nReaches, 1);
for iReach = 1:nReaches
    % get (filtered) accel data around each reach
    isPeriReach = (tst > (reachOnset(iReach))) & ...
        (tst <= (reachEnd(iReach) + 0.5));
    accReach = accMag(isPeriReach,:);
    t = (1/fs) * (0:(length(accReach)-1));
    
    % Get summed instantaneous triaxial power in that band
%     accMag = abs(hilbert(accReach));
    
%     figure; 
%     plot(t, accReach);
%     figure; 
%     plot(t, accMag);
%     figure; 
%     plot(t, sum(accReach, 2))
    reachPow(iReach) = mean(sum(accReach, 2));

    
end

figure; violinplot(reachPow);






