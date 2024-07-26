% script to convert from TDTblock to NWB file
% 
% 
% % TDTblock = 'AzNeuralverTwo-210726-105657';
% TDTblock = 'AzNeuralverTwo-210726-105928';
% % TDTblock = 'AzNeuralverTwo-210726-110153';
% % TDTblock = 'AzNeuralverTwo-210726-110428';
% % TDTblock = 'AzNeuralverTwo-210726-110705';
% % TDTblock = 'AzNeuralverTwo-210726-110934';
% % TDTblock = 'AzNeuralverTwo-210726-111240';
% % TDTblock = 'AzNeuralverTwo-210726-111525';
% % TDTblock = 'AzNeuralverTwo-210726-111915';
% % TDTblock = 'AzNeuralverTwo-210726-112157';
% % TDTblock = 'AzNeuralverTwo-210726-112448';
% % TDTblock = 'AzNeuralverTwo-210726-112727';
% 
% 
% 
% 
% % Read in basic stream info from TDT
% pn = 'C:\DATAtemp\STNstim_GPrecord\Data Acquisition\Nebula-210726\';
% % TDTblock = 'AzNeuralverTwo-210726-105657';
% 
% blockpath = [pn TDTblock];
% 
% 
% nwb = contrib.tdt.TDTbin2NWB([pn TDTblock])
% 
% % write and read back an NWB file
% % write file
% nwbfn = nwb.general_session_id;
% nwbfullpath = [pn nwbfn '.nwb'];
% disp(['Exporting ' nwbfn '.nwb...']);
% tic
% nwbExport(nwb, nwbfullpath)
% toc
% disp('DONE!')
% 
% 
% % test read
% nwbTest = nwbRead(nwbfullpath, 'ignoreCache');


%% Test the case where I specify the data conversion of each TDT stream


% script to convert from TDTblock to NWB file


% TDTblock = 'AzNeuralverTwo-210726-105657';
% TDTblock = 'AzNeuralverTwo-210726-105928';
% TDTblock = 'AzNeuralverTwo-210726-110153';
TDTblock = 'AzNeuralverTwo-210726-110428';
% TDTblock = 'AzNeuralverTwo-210726-110705';
% TDTblock = 'AzNeuralverTwo-210726-110934';
% TDTblock = 'AzNeuralverTwo-210726-111240';
% TDTblock = 'AzNeuralverTwo-210726-111525';
% TDTblock = 'AzNeuralverTwo-210726-111915';
% TDTblock = 'AzNeuralverTwo-210726-112157';
% TDTblock = 'AzNeuralverTwo-210726-112448';
% TDTblock = 'AzNeuralverTwo-210726-112727';




% Read in basic stream info from TDT
pn = 'C:\DATAtemp\STNstim_GPrecord\Data Acquisition\Nebula-210726\';
% TDTblock = 'AzNeuralverTwo-210726-110153';

blockpath = [pn TDTblock];

streamParams(1).name = 'Wav2';
streamParams(1).data_unit = 'Volts';
streamParams(1).data_conversion = 1/(1e9);
streamParams(2).name = 'pSUg';
streamParams(2).data_unit = '<missing>';
streamParams(2).data_conversion = 1;
streamParams(3).name = 'Dbs1';
streamParams(3).data_unit = '<missing>';
streamParams(3).data_conversion = 1;
streamParams(4).name = 'RSn1';
streamParams(4).data_unit = '<missing>';
streamParams(4).data_conversion = 1;

nwb = contrib.tdt.TDTbin2NWB([pn TDTblock], streamParams);

% write and read back an NWB file
% write file
nwbfn = nwb.general_session_id;
nwbfullpath = [pn nwbfn '.nwb'];
disp(['Exporting ' nwbfn '.nwb...']);
tic
nwbExport(nwb, nwbfullpath)
toc
disp('DONE!')


% test read
nwbTest = nwbRead(nwbfullpath, 'ignoreCache');