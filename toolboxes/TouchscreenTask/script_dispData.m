% script for reading in the excel sheet with TS task details, then
% displaying charts



% Read in excel sheet to Table
pn = 'K:\My Drive\PROJECTS\TouchscreenTask\TrainingRecords';
fn = 'TouchScreenTrainingRecord_Uva.xlsx';

T = readtable([pn, '\', fn]);



% Convert values to datetime arrays



% plot average reaches per week! 