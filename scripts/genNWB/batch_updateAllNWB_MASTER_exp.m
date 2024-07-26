% batch script for updating all NWB master files (main reference file; other
% NWB data not contained within the master file has reference "external
% links" within the master file, so everything can be called/accessed by
% interfacing with the master file. 

% When running this updated script for the master files, code assumes that:
% 1) All new data to be included as extra NWB files has already been
% generated and that external links are already figured out within the code
% "script_genNWB_MASTER_exp.m". 
% 2) The latest versions of teh NWB master files are either deleted or
% moved away from the target directory in which the current code will save
% the newest updated NWB master files. NWB functions don't overwrite old
% files of the same name by default, that situation makes the code break,
% at least in the current version that Ben Dichter and co. have written...

% run the NWB master creation script for each harmaline experiment day:

clear;
SESSION_ID = 'TremoLfpDBS-190927-100155';
script_genNWB_MASTER_exp

clear;
SESSION_ID = 'TremoLfpDBS-191004-100637';
script_genNWB_MASTER_exp

clear;
SESSION_ID = 'TremoLfpDBS-191011-104322';
script_genNWB_MASTER_exp

clear;
SESSION_ID = 'TremoLfpDBS-191018-100615';
script_genNWB_MASTER_exp

clear;
SESSION_ID = 'TremoLfpDBS-191025-104651';
script_genNWB_MASTER_exp

clear;
SESSION_ID = 'TremoLfpDBS-191101-101430';
script_genNWB_MASTER_exp

clear;
SESSION_ID = 'TremoLfpDBS-191108-101829';
script_genNWB_MASTER_exp

clear;
SESSION_ID = 'TremoLfpDBS-191115-100127';
script_genNWB_MASTER_exp