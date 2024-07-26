% script to plot out the coefficients of the formula fit to washin data





% load in tabulated values, and specify the coefficient of interst
projRootPn = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\';
coeffTabPn = 'Data Processing\washinFits\';
T = readtable([projRootPn coeffTabPn 'betaLowFitCoeffs.xlsx']);

coeffStr = 'k';

% re-refernce data for errorbar plotting
coeffData = [T{:,coeffStr}, T{:,[coeffStr 'lower']}, T{:,[coeffStr 'upper']}];
coeffReref = coeffData;
% re-reference "lower" to be the length below the average
coeffReref(:,2) = coeffData(:,1) - coeffData(:,2);

% re-reference "uppper" to be the length above the average
coeffReref(:,3) = coeffData(:,3) - coeffData(:,1);

% display results
figure; ax = axes;
e = errorbar([1:8], coeffReref(:,1), coeffReref(:,2), coeffReref(:,3)); 
e.LineStyle = 'none';
e.Marker = 'o';
labels = {''; '1st 2 mg/kg'; '1st 4 mg/kg'; '1st 6 mg/kg'; '1st 8 mg/kg'; '2nd 2 mg/kg'; '2nd 4 mg/kg'; '2nd 6 mg/kg'; '2nd 8 mg/kg'; ''; ''}; 
ax.XLim = [0, 10];
% ax.YLim = [2.5, 7];
ax.XTickLabel = labels;
ax.XTickLabelRotation = -45;

ylabel(['coefficient:   ' coeffStr]);
title(['Exponential function fits to Washin Low Beta'])