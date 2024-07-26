function nwb = func_removeArtifacts_NWB(nwb, cfg)
% detect DBS stim times based on the "Dbs1" store in the TDT
% block recorded, now in nwb file format. Then perform artifact subtraction
% 
% INPUTS
% 
% nwb:
% NWBfile object
%
% cfg. 
%    .threshold: (scalar)
%       point at which DBS aftifacts cross
%
%    .procSpkChans: (1xN vector)
%       1xN vector of channel subselection to perform processing on and 
%       save under /processing in nwb file.
%
%    .artDetect_offest: (scalar)
%       number of samples to shift (right) the detected artifact onset; 
%       since artifacts were detected on Dbs1, but subtraction applied on 
%       Wav2, this correction helps with the sample-offset problem.
%
%    .artDetect_blankWindow: (1x2 vector)
%       left and right edges of a window to set data to zero around each 
%       DBS-pulse; number of samples.
%
%    .spkBandPass: (1x2 vector)
%       detailing bandpass frequency limits for a butterworth filter.
%
% OUTPUTS
% 
% nwbUpd:
% NWBfile object with all processed data added (as TimeSeries objects)
% under /procssing within the "DBSartifact_detect-clean" ProcessingModule
% object. 
% spkClean: Artifact-cleaned and spike-filtered continuous data for
% specified cfg.procSpkChans.
% tDBS: times in seconds of detected DBS artifacts
% tDBS_uncorrected: same as tDBS, but before the offset-correction step.
% tDBS_virtPre: artificially-generated "virtual" DBS pulses occurring in
% all data prior to the onset of the first real detected DBS pulse. 
% tDBS_virtPos: same as above, but for period after final detected DBS
% pulse.



Dbs1 = nwb.acquisition.get('Dbs1');
dbsTrace = double(Dbs1.data.load);
fs_res = Dbs1.starting_time_rate;
tStamp = (1/fs_res) * (0:(length(dbsTrace) - 1));
% figure; ax = axes; plot(ax, tStamp, dbsTrace)



%% Detect DBS times based on the voltage trace Dbs1

% THRESH = -0.01; % manually change this basded on observation
THRESH = cfg.threshold; % manually change this basded on observation

x = dbsTrace; 

if THRESH < 0 % if thresh is a negative crossing, flip the data over 
    x = -x; 
    THRESH = -THRESH;
    
end

% Essentially need to pick out values which exceed threshold with the condition that the previous value  
% needs to be below the threshold
idxl = x >= THRESH;
idxl(1) = 0;
idx = find(idxl);
yest = x(idx-1) < THRESH; 
idxStim = idx(yest); % Final output
idxStim = idxStim - 1; % correct by shifting stim onset one sample back
stimTime_uncorr = tStamp(idxStim);

disp(['detected stims: ' num2str(length(stimTime_uncorr))]);
% should be just one bin; more bins means that multiple parts of the pulse
% are getting detected 
figure; histogram(diff(stimTime_uncorr));



%% Load relevant raw data contaminated with DBS pulse artifacts 
% Specify the subset of recorded channels that are relevant to spike
% analysis for this recording session. Upsample them if needed to match the
% sampling frequency of the Dbs1 channel.

% chans = [64, 65, 66, 68, 70, 71, 74, 78];
chans = cfg.procSpkChans;



% sampFactor = 2;
tSer = nwb.acquisition.get(cfg.spkStore);

% Get original timestamps
N = tSer.data.dims(1);
t0 = tSer.starting_time;
fs = tSer.starting_time_rate;
tst = t0 + ((0:(N-1))) * (1/fs);

% % Get Resampled timestamps
% N_res = N * sampFactor;
% t0_res = t0 / sampFactor;
% fs_res = fs * sampFactor;
% tst_res = ((1/fs_res) * (0:(N_res-1))) + t0_res;
% 
% nChans = length(chans);
% raw = zeros(N_res, nChans);
% for iCh = 1:nChans
%     data = double(tSer.data.load_mat_style(:,chans(iCh)));
%     data_res = interp1(tst, data, tst_res);
%     data_res(1) = data_res(2);
%     data_res(end) = data_res(end-1); % correct for any NaNs output by interp1...
%     raw(:,iCh) = data_res;
% 
% end
% 
nChans = length(chans);
raw = zeros(N, nChans);
for iCh = 1:nChans
    data = double(tSer.data.load_mat_style(:,chans(iCh)));
%     data_res = interp1(tst, data, tst_res);
%     data_res(1) = data_res(2);
%     data_res(end) = data_res(end-1); % correct for any NaNs output by interp1...
%     raw(:,iCh) = data_res;
    raw(:,iCh) = data;


end



%% Remove DBS artifact using a subtraction method

% CONSTANTS

% Peri-stimuls windowing and correction
% OFFSET_CORRECT = 27; % samples; if detected stim times were several samples to the left or right of where the pulse appears in data, correct it here...
OFFSET_CORRECT = cfg.artDetect_offest; % samples; if detected stim times were several samples to the left or right of where the pulse appears in data, correct it here...

% blankWin = [-4, 8]; % sample-window for blanking around stim pulse detection...
blankWin = cfg.artDetect_blankWindow; % sample-window for blanking around stim pulse detection...


% Spike-filtering
% FC = [200, 8000]; % Hz, bandpass for butterworth filter...
FC = cfg.spkBandPass; % Hz, bandpass for butterworth filter...

stimSamps_uncorr = round(stimTime_uncorr * fs_res);
stimSamps = stimSamps_uncorr + OFFSET_CORRECT; % offset by 4 samps to the left, as was done in SARGE. 
stimTime = stimSamps / fs_res;
% loffset = 0;
% nTempl = 30;

params.fc = FC; % Hz, bandpass filter

% Temporarily shift stim-detection times and total blanked samples so
% taht the desired peri-stim blank window is achieved, then remove
% Artifacts:
params.blankSamps = blankWin(2) - blankWin(1); % n samples to right or stim samples to blank

spk = raw;
for iCh = 1:nChans
    spk(:,iCh) = remArt(raw(:,iCh), fs_res, stimSamps + blankWin(1), params);
        
end



%% Create virtual pulse time events for the pre- and post- DBS portion of data

isi = median(diff(stimSamps));

% post-DBS virt pulse times
stim0 = stimSamps(end) + isi;
virtPos_idx = stim0:isi:size(spk, 1);
virtPos_idx(end) = []; % remove last index, so theres froom for blanking
virtPos_tst = tst(virtPos_idx);

% pre-DBS virt pulse times
stimEnd = stimSamps(1) - isi;
virtPre_idx = 1:isi:stimEnd;
virtPre_idx(1) = []; % remove first index, so theres room for blanking
virtPre_tst = tst(virtPre_idx);

% Create zero-blanked regions around all pulses, virtual or real
blankMask = ones(size(spk, 1), 1);

for i = 1:length(virtPre_idx)
    b0 = virtPre_idx(i) + blankWin(1);
    bend = virtPre_idx(i) + blankWin(2);
    blankMask(b0:bend) = 0;
    
end

for i = 1:length(stimSamps)
    b0 = stimSamps(i) + blankWin(1);
    bend = stimSamps(i) + blankWin(2);
    blankMask(b0:bend) = 0;
    
end

for i = 1:length(virtPos_idx)
    b0 = virtPos_idx(i) + blankWin(1);
    bend = virtPos_idx(i) + blankWin(2);
    blankMask(b0:bend) = 0;
    
end

blankMask2 = repmat(blankMask, 1, nChans);

spkClean = spk .* blankMask2;



%% Save updates in the NWB file
figure; plot(tStamp, dbsTrace, tStamp, spkClean);
% Update the current nwb file to include DBS stim times and clean data

tDBS_uncorrected = types.core.TimeSeries(...
    'data', ones(length(stimTime_uncorr), 1), ...
    'data_unit', 'N/A', ...
    'data_continuity', 'instantaneous', ...
    'description', 'DBS event timings found in timestamps field, the data field is just placeholder values', ...
    'timestamps', stimTime_uncorr);

tDBS_virtPre = types.core.TimeSeries(...
    'data', ones(length(virtPre_tst), 1), ...
    'data_unit', 'N/A', ...
    'data_continuity', 'instantaneous', ...
    'description', 'DBS event timings found in timestamps field, the data field is just placeholder values', ...
    'timestamps', virtPre_tst);

tDBS = types.core.TimeSeries(...
    'data', ones(length(stimTime), 1), ...
    'data_unit', 'N/A', ...
    'data_continuity', 'instantaneous', ...
    'description', 'DBS event timings found in timestamps field, the data field is just placeholder values', ...
    'timestamps', stimTime);

tDBS_virtPos = types.core.TimeSeries(...
    'data', ones(length(virtPos_tst), 1), ...
    'data_unit', 'N/A', ...
    'data_continuity', 'instantaneous', ...
    'description', 'DBS event timings found in timestamps field, the data field is just placeholder values', ...
    'timestamps', virtPos_tst);


TSspkSclean = types.core.TimeSeries( ...
        'data', spkClean, ...
        'data_continuity', tSer.data_continuity, ...
        'data_conversion', tSer.data_conversion, ...
        'data_resolution', tSer.data_resolution, ...
        'data_unit', tSer.data_unit, ...
        'description', 'Clean data for analysis includes chans: 64, 65, 66, 68, 70, 71, 74, 78', ...
        'starting_time', t0, ... % seconds
        'starting_time_rate', fs);


ecephys_module = types.core.ProcessingModule(...
    'description', 'Holds ephys data that has been DBS artifact-cleaned');
ecephys_module.nwbdatainterface.set('spkClean', TSspkSclean);
ecephys_module.nwbdatainterface.set('tDBS', tDBS);
ecephys_module.nwbdatainterface.set('tDBS_virtPre', tDBS_virtPre);
ecephys_module.nwbdatainterface.set('tDBS_virtPos', tDBS_virtPos);
ecephys_module.nwbdatainterface.set('tDBS_uncorrected', tDBS_uncorrected);

nwb.processing.set('DBSartifact_detect-clean', ecephys_module);



end

