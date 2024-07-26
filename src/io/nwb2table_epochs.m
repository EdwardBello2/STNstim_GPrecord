function tbl = nwb2table_epochs(nwb)
% function for converting the nwb dynamic table to matlab table, specific
% to the interval-epochs information from the nwb file. Time data is in
% seconds, relative to session start time in the nwb file by definition. 

dyntbl = nwb.intervals_epochs;
id = dyntbl.id.data.load;
varDescr{1} = 'record id';

tbl = table(id);

colnames = dyntbl.colnames;
nColnames = numel(colnames);
for iCol = 1:nColnames
    eval([colnames{iCol} ' = dyntbl.' colnames{iCol} '.data.load;']);
    eval(['tbl = [tbl, table(' colnames{iCol} ')];']);
      
    varDescr{iCol+1} = dyntbl.(colnames{iCol}).description;
    
end


tbl = sortrows(tbl, 'start_time');


tbl.Properties.Description = dyntbl.description;
% tbl.Properites.VariableDescriptions = ;
% tbl.Properites.VariableUnits = {''; ''; 'seconds'; 'seconds'};



end