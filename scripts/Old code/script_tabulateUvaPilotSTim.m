


%% Load metadata for all files

% Read in metadata table for Mela's recordings
projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';
tablePath = 'Data Acquisition\Uva Pilot\';
metaTab = readtable([projRootPath tablePath 'acquisitionMetadata_Uva.csv']);

dataAcqPn = 'D:\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TDTdata\';
% directoryPn = 'Mela\Harmaline\';


% harTab = readtable([projRootPath tablePath 'harmalineSessionData.csv']);
% fullMeta = join(metaTab, harTab);
% 
% % add in data for time since harmaline injection
% harRefTime = fullMeta.startTime - fullMeta.harInjTime;
% fullMeta = [fullMeta, table(harRefTime)];
% 
% % Keep only rows (recordings) marked as "ok" in recComment
% fullMeta = fullMeta(strcmp(fullMeta.recComment, 'ok'),:);

% Create table to track DBS stim parameters; particularly which channel was
% disconnected from recording so that DBS could be delivered thru it

stimTab = metaTab(:,1:2); % Base table off of metaTab



% load epocs field for each tdt block and collect the stim parameter data
% for taht file

nRows = height(stimTab);

     period = zeros(nRows, 1);
        amp = zeros(nRows, 1);
      tdtCh = zeros(nRows, 1);
elecContact = cell(nRows, 1);
     nStims = zeros(nRows, 1);
  timeOnset = zeros(nRows, 1);
  
  elecLabels = {'C0'; 'C1'; 'C2'; 'C3'; 'C4'; 'C5'; 'C6'; 'C7'};

for iRow = 1:nRows
    % load epocs
    tdt = TDTbin2mat([ dataAcqPn stimTab.TDTtank{iRow} '\' stimTab.TDTblock{iRow}], ...
        'TYPE', {'epocs'});
    
    
    if isfield(tdt.epocs, 'Pper')
        if numel(tdt.epocs.Pper.data) > 1, error('stim blocks > 1'); end
        period(iRow,1) = tdt.epocs.Pper.data;
        
    else
        period(iRow,1) = NaN;
        
    end
    
    
    if isfield(tdt.epocs, 'AmA_')
        if ~all(tdt.epocs.AmA_.data == tdt.epocs.AmA_.data(1))
            error('stim blocks > 1');
            
        end
        amp(iRow,1) = tdt.epocs.AmA_.data(1);
        
    else
        amp(iRow,1) = NaN;
        
    end
        
    
    if isfield(tdt.epocs, 'StCh')
        if ~all(tdt.epocs.StCh.data == tdt.epocs.StCh.data(1))
            error('stim blocks > 1');
            
        end
        tdtCh(iRow,1) = tdt.epocs.StCh.data(1);
        elecContact{iRow,1} = elecLabels{tdt.epocs.StCh.data(1)};
        
    else
        tdtCh(iRow,1) = NaN;
        elecContact{iRow,1} = 'missing';
        
    end
          
            
    if isfield(tdt.epocs, 'Pct_')
        if numel(tdt.epocs.Pct_.data) > 1, error('stim blocks > 1'); end
        nStims(iRow,1) = tdt.epocs.Pct_.data;
        
    else
        nStims(iRow,1) = NaN;
        
    end
    
    if isfield(tdt.epocs, 'Pper')
        if numel(tdt.epocs.Pper.onset) > 1, error('stim blocks > 1'); end
        timeOnset(iRow,1) = tdt.epocs.Pper.onset;
        
    else
        timeOnset(iRow,1) = NaN;
        
    end

end
 
  
stimTab = [stimTab, table(period), table(amp), table(tdtCh), ...
    table(elecContact), table(nStims), table(timeOnset)];

writetable(stimTab, [projRootPath tablePath 'dbsMetaData_Uva.csv']);


%% SUB-FUNCTIONS





