function nwb = nwbRead_specSess(cfg, sessionMetaData)
% load the matlab nwb object from one of the nwb files from the harmaline
% experiment. Calls nwbRead.m from the matnwb toolbox
% isTest = false;

pn = cfg.projRootDirFullpath;
pn = assertEndSlash(pn);

nwbPath = 'Data Processing\NWB\';

% Get basic metadata table
T = readtable([pn nwbPath 'nwbfile_basicMetadata.xlsx']);


cols.harmalineDose = sessionMetaData.dose;
cols.exposure = sessionMetaData.exposureNum;
cols.sessionType = sessionMetaData.sessionType;

subT = getRows(T, cols);
nwbfn = subT.nwbFile_name{1};
% if isTest, nwbfn = [nwbfn 'TEST']; end
    


nwb = nwbRead([pn nwbPath nwbfn], 'ignorecache'); % ignorecache makes sure a +types folder doesn't get made again in current folder

if isempty(nwb.general_session_id), nwb.general_session_id = nwbfn; end

end

%% SUB-FUNCTIONS

function subT = getRows(T, cols)

colNames = fieldnames(cols);
nFields = length(fieldnames(cols));

isIdx = false(height(T), nFields);
for iF = 1:nFields
    isIdx(:,iF) = strcmp(cols.(colNames{iF}), T{:,colNames{iF}});
    
end

select = all(isIdx, 2);
subT = T(select,:);


end