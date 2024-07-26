% test for reading in data from BasicParadigmEnhanced touchscreen task
% software, number of successful reaches being focus

% All data comes from the spreadsheet "metadataByDay.xlsx" in 
% "Thalamic DBS for Harmaline Tremors\Data Acquisition\TStask BasicParadigmEnhanced\harmalineDoseStudy"

% Display grouped barplots
% data = {rand(100,4), rand(20,4)*.8, rand(1000,4)*1.2};
% figure; axes;
% boxplotGroup(data)

% har = [53 121; 26 13; 15 92; 0 1;]; % Fill manually with successTrial counts from each harmaline day
% m = []; % fill this manually with just maintenance day successTrial data

%% load in TS taks session metadata
clear; 

tabFn = 'metadataByDay_NAremoved.xlsx';
pn = 'L:\My Drive\PROJECTS\Thalamic DBS for Harmaline Tremors\Data Acquisition\TStask BasicParadigmEnhanced\harmalineDoseStudy\';
T = readtable([pn tabFn]);



%% Calc and add individual trial reaction-times and reach-times

% Loop thru all days
nRows = height(T);
reactionTime = cell(nRows, 1);
reachTime = cell(nRows, 1);
reactionTime_av = zeros(nRows, 1);
reachTime_av = zeros(nRows, 1);

for iRow = 1:nRows
    % Read in text from a given day
    fn = T.filename{iRow};
    inChar = fileread([pn fn '\' fn '.txt']);
    [~, ReachTrials] = parseInfoRound(inChar);

    % Get all reaction-times and reach times for each day
    reactionTime{iRow,1} = ReachTrials.ReactionTime;
    reactionTime_av(iRow,1) = mean(ReachTrials.ReactionTime);
    
    reachTime{iRow,1} = ReachTrials.ReachTime;
    reachTime_av(iRow,1) = mean(ReachTrials.ReachTime);

end
T = [T, table(reactionTime), table(reactionTime_av), ...
    table(reachTime), table(reachTime_av)];



%% Show distributions of session reaction-times 

% Harmaline-session Reaction-times
rxnT.exp1dos2 = T.reactionTime{(T.harDose == 2) & (T.doseExpos == 1) & (T.harDelivered == true)};
rxnT.exp2dos2 = T.reactionTime{(T.harDose == 2) & (T.doseExpos == 2) & (T.harDelivered == true)};
rxnT.exp1dos4 = T.reactionTime{(T.harDose == 4) & (T.doseExpos == 1) & (T.harDelivered == true)};
rxnT.exp2dos4 = T.reactionTime{(T.harDose == 4) & (T.doseExpos == 2) & (T.harDelivered == true)};
rxnT.exp1dos6 = T.reactionTime{(T.harDose == 6) & (T.doseExpos == 1) & (T.harDelivered == true)};
rxnT.exp2dos6 = T.reactionTime{(T.harDose == 6) & (T.doseExpos == 2) & (T.harDelivered == true)};
% rxnT.exp1dos8 = T.reactionTime{(T.harDose == 8) & (T.doseExpos == 1) & (T.harDelivered == true)};
% rxnT.exp2dos8 = T.reactionTime{(T.harDose == 8) & (T.doseExpos == 2) & (T.harDelivered == true)};

% % remove outlier data from all harmaline days
% rxnT.exp1dos2 = util.remoutliers(rxnT.exp1dos2, 'MADthresh', 3);
% rxnT.exp2dos2 = util.remoutliers(rxnT.exp2dos2, 'MADthresh', 3);
% rxnT.exp1dos4 = util.remoutliers(rxnT.exp1dos4, 'MADthresh', 3);
% rxnT.exp2dos4 = util.remoutliers(rxnT.exp2dos4, 'MADthresh', 3);
% rxnT.exp1dos6 = util.remoutliers(rxnT.exp1dos6, 'MADthresh', 3);
% rxnT.exp2dos6 = util.remoutliers(rxnT.exp2dos6, 'MADthresh', 3);



% naive-session Reaction-times
fieldStr = 'reactionTime';

% Get maintenance day (naive) data, and remove outlier values
mOutliers = cell2mat(T.(fieldStr)(T.harDelivered == false)); 
m = util.remoutliers(mOutliers, 'MADthresh', 3);
ci = prctile(m, 100*[0.05, 0.95]);

% m = cell2mat(T.(fieldStr)(T.harDelivered == false)); % fill this manually with just maintenance day successTrial data


figure;
% Display violin plot of maintenance days
ax(1) = subplot(1,6,1); 
violinplot(m, [], 'Width', 0.3, 'ShowData', false);
grid on; hold on
plot(ax(1), [ax(1).XLim], [ci(1) ci(1)], '--r', [ax(1).XLim], [ci(2) ci(2)], '--r');
title('Mean reaction-time of trials on all Naive-maintenace TS task days'); 
% ylabel('# target-reaches'); grid on;

% Display grouped violin plots of harmaline days
ax(2) = subplot(1,6,2:6);
violinplot(rxnT); 
title('Individual trial reaction-times on specific Harmaline exposure TS task days');
% ax(2).XTickLabel = {'2 mg/kg'; '4 mg/kg'; '6 mg/kg'; '8 mg/kg'};
grid on; hold on
plot(ax(2), [ax(2).XLim], [ci(1) ci(1)], '--r', [ax(2).XLim], [ci(2) ci(2)], '--r');


set(ax(:), 'YLim', [0 3]);
set(gcf, 'Position', [2107 546 1371 420]);



%% Show distributions of session reach-times

% Harmaline-session Reach-times
rchT.exp1dos2 = T.reachTime{(T.harDose == 2) & (T.doseExpos == 1) & (T.harDelivered == true)};
rchT.exp2dos2 = T.reachTime{(T.harDose == 2) & (T.doseExpos == 2) & (T.harDelivered == true)};
rchT.exp1dos4 = T.reachTime{(T.harDose == 4) & (T.doseExpos == 1) & (T.harDelivered == true)};
rchT.exp2dos4 = T.reachTime{(T.harDose == 4) & (T.doseExpos == 2) & (T.harDelivered == true)};
rchT.exp1dos6 = T.reachTime{(T.harDose == 6) & (T.doseExpos == 1) & (T.harDelivered == true)};
rchT.exp2dos6 = T.reachTime{(T.harDose == 6) & (T.doseExpos == 2) & (T.harDelivered == true)};
% rchT.exp1dos8 = T.reachTime{(T.harDose == 8) & (T.doseExpos == 1) & (T.harDelivered == true)};
% rchT.exp2dos8 = T.reachTime{(T.harDose == 8) & (T.doseExpos == 2) & (T.harDelivered == true)};

% naive-session reach-times
fieldStr = 'reachTime';

% Get maintenance day (naive) data, and remove outlier values
mOutliers = cell2mat(T.(fieldStr)(T.harDelivered == false)); 
m = util.remoutliers(mOutliers, 'MADthresh', 3);
ci = prctile(m, 100*[0.05, 0.95])

% m = cell2mat(T.(fieldStr)(T.harDelivered == false)); % fill this manually with just maintenance day successTrial data


figure;
% Display violin plot of maintenance days
ax(1) = subplot(1,6,1); 
violinplot(m, [], 'Width', 0.3, 'ShowData', false);
grid on; hold on
plot(ax(1), [ax(1).XLim], [ci(1) ci(1)], '--r', [ax(1).XLim], [ci(2) ci(2)], '--r');
title('Mean reach-time of trials on all Naive-maintenace TS task days'); 
% ylabel('# target-reaches'); grid on;

% Display grouped violin plots of harmaline days
ax(2) = subplot(1,6,2:6);
violinplot(rchT);
title('Individual trial reach-times on specific Harmaline exposure TS task days');
% ax(2).XTickLabel = {'2 mg/kg'; '4 mg/kg'; '6 mg/kg'; '8 mg/kg'};
grid on; hold on;
plot(ax(2), [ax(2).XLim], [ci(1) ci(1)], '--r', [ax(2).XLim], [ci(2) ci(2)], '--r');


set(ax(:), 'YLim', [0 3]);
set(gcf, 'Position', [2106 42 1371 420]);



%% Show all SuccessTrials

fieldStr = 'SuccessTrials';

har = [T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 1) & (T.harDelivered == true)), ... 
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 2) & (T.harDelivered == true)); ... 
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 2) & (T.harDelivered == true));];

m = T.(fieldStr)(T.harDelivered == false); % fill this manually with just maintenance day successTrial data

figure;
% Display violin plot of maintenance days
ax(1) = subplot(1,3,1); 
violinplot(m, [], 'Width', 0.3);
title('Naive-maintenace TS task days'); 
ylabel('# target-reaches'); grid on;

% Display grouped barplots of harmaline days
ax(2) = subplot(1,3,2:3);
bar(har); legend('expos 1', 'expos 2');
title('Harmaline exposure TS task days');
ax(2).XTickLabel = {'2 mg/kg'; '4 mg/kg'; '6 mg/kg'; '8 mg/kg'};
grid on;

set(ax(:), 'YLim', [0 220]);
set(gcf, 'Position', [2008 360 842 420]);



%% Show all ReplicationsCompleted

fieldStr = 'ReplicationsCompleted';

har = [T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 1) & (T.harDelivered == true)), ... 
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 2) & (T.harDelivered == true)); ... 
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 2) & (T.harDelivered == true));];

m = T.(fieldStr)(T.harDelivered == false); % fill this manually with just maintenance day successTrial data

figure;
% Display violin plot of maintenance days
ax(1) = subplot(1,3,1); 
violinplot(m, [], 'Width', 0.3);
title('Naive-maintenace TS task days'); 
grid on
% ylabel('# target-reaches'); grid on;

% Display grouped barplots of harmaline days
ax(2) = subplot(1,3,2:3);
bar(har); legend('expos 1', 'expos 2');
title('Harmaline exposure TS task days');
ax(2).XTickLabel = {'2 mg/kg'; '4 mg/kg'; '6 mg/kg'; '8 mg/kg'};
grid on;

set(ax(:), 'YLim', [0 240]);
set(gcf, 'Position', [2008 360 842 420]);



%% Show number of errors-reaches

fieldStr = 'ErrorTrials';

har = [T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 1) & (T.harDelivered == true)), ... 
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 2) & (T.harDelivered == true)); ... 
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 2) & (T.harDelivered == true));];

m = T.(fieldStr)(T.harDelivered == false); % fill this manually with just maintenance day successTrial data

figure;
% Display violin plot of maintenance days
ax(1) = subplot(1,3,1); 
violinplot(m, [], 'Width', 0.3);
title('Naive-maintenace TS task days'); 
% ylabel('# target-reaches'); grid on;

% Display grouped barplots of harmaline days
ax(2) = subplot(1,3,2:3);
bar(har); legend('expos 1', 'expos 2');
title('Harmaline exposure TS task days');
ax(2).XTickLabel = {'2 mg/kg'; '4 mg/kg'; '6 mg/kg'; '8 mg/kg'};
grid on;

set(ax(:), 'YLim', [0 20]);
set(gcf, 'Position', [2008 360 842 420]);



%% Show percentage of error reaches relative to all replications

fieldStr = 'x_errorReach';

har = [T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 1) & (T.harDelivered == true)), ... 
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 2) & (T.harDelivered == true)); ... 
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 2) & (T.harDelivered == true));];

m = T.(fieldStr)(T.harDelivered == false); % fill this manually with just maintenance day successTrial data

figure;
% Display violin plot of maintenance days
ax(1) = subplot(1,3,1); 
violinplot(m, [], 'Width', 0.3);
title('Naive-maintenace TS task days'); 
% ylabel('# target-reaches'); grid on;

% Display grouped barplots of harmaline days
ax(2) = subplot(1,3,2:3);
bar(har); legend('expos 1', 'expos 2');
title('Harmaline exposure TS task days');
ax(2).XTickLabel = {'2 mg/kg'; '4 mg/kg'; '6 mg/kg'; '8 mg/kg'};
grid on;

set(ax(:), 'YLim', [0 70]);
set(gcf, 'Position', [2008 360 842 420]);



%% Show number of Premature Reach-initiations

fieldStr = 'InvalidInit';

har = [T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 1) & (T.harDelivered == true)), ... 
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 2) & (T.harDelivered == true)); ... 
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 2) & (T.harDelivered == true));];

m = T.(fieldStr)(T.harDelivered == false); % fill this manually with just maintenance day successTrial data

figure;
% Display violin plot of maintenance days
ax(1) = subplot(1,3,1); 
violinplot(m, [], 'Width', 0.3);
title('Naive-maintenace TS task days'); 
% ylabel('# target-reaches'); grid on;

% Display grouped barplots of harmaline days
ax(2) = subplot(1,3,2:3);
bar(har); legend('expos 1', 'expos 2');
title('Harmaline exposure TS task days');
ax(2).XTickLabel = {'2 mg/kg'; '4 mg/kg'; '6 mg/kg'; '8 mg/kg'};
grid on;

set(ax(:), 'YLim', [0 125]);
set(gcf, 'Position', [2008 360 842 420]);



%% Total number of Reach attempts (Premature Reaches + Proper Reaches of any kind that connected)

fieldStr = 'ReachAttempts';

har = [T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 1) & (T.harDelivered == true)), ... 
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 2) & (T.harDelivered == true)); ... 
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 2) & (T.harDelivered == true));];

% Get maintenance day (naive) data, and remove outlier values
mOutliers = T.(fieldStr)(T.harDelivered == false); 
m = util.remoutliers(mOutliers, 'MADthresh', 3);
ci = prctile(m, 100*[0.05, 0.95])


figure;
set(gcf, 'Position', [2008 360 842 420]);

Ymax = 360;

% Display violin plot of maintenance days
ax(1) = subplot(1,3,1); 
violinplot(m, [], 'Width', 0.3);
grid on; hold on;
title('Naive-maintenace TS task days'); 
plot(ax(1), [ax(1).XLim], [ci(1) ci(1)], '--r', [ax(1).XLim], [ci(2) ci(2)], '--r');

% ylabel('# target-reaches'); grid on;

% Display grouped barplots of harmaline days
ax(2) = subplot(1,3,2:3);
set(ax(2), 'XLim', [0.5143 4.4857]);
plot(ax(2), [ax(2).XLim], [ci(1) ci(1)], '--r', [ax(2).XLim], [ci(2) ci(2)], '--r');
hold on
bar(har); legend(ax(2).Children(1:2), 'expos 2', 'expos 1');
title('Harmaline exposure TS task days');
ax(2).XTickLabel = {'2 mg/kg'; ''; '4 mg/kg'; ''; '6 mg/kg'; ''; '8 mg/kg'};
grid on;

set(ax(:), 'YLim', [0 Ymax]);


% Get the 95% percentile intervals of the data (assuming a normal
% distribution)

% ci = norminv([0.025, 0.975], mean(m),  std(m));
% bootfun = @(x)(x);
% ci = prctile(m, 100*[0.05, 0.95])
% plot(ax(1), [ax(1).XLim], [ci(1) ci(1)], '--r', [ax(1).XLim], [ci(2) ci(2)], '--r');
% plot(ax(2), [ax(2).XLim], [ci(1) ci(1)], '--r', [ax(2).XLim], [ci(2) ci(2)], '--r');






%% % of all Reach attempts that were premature-reaches

fieldStr = 'x_invalidInitVSallReachAttempts';

har = [T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 1) & (T.harDelivered == true)), ... 
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 2) & (T.harDelivered == true)); ... 
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 2) & (T.harDelivered == true));];

% Get maintenance day (naive) data, and remove outlier values
mOutliers = T.(fieldStr)(T.harDelivered == false); 
m = util.remoutliers(mOutliers, 'MADthresh', 3);
ci = prctile(m, 100*[0.05, 0.95]); % 90% interval


figure;
set(gcf, 'Position', [2007 90 986 420]);

Ymax = 100;

% Display violin plot of maintenance days
ax(1) = subplot(1,3,1); 
violinplot(m, [], 'Width', 0.3);
grid on; hold on;
title('Naive-maintenace TS task days'); 
plot(ax(1), [ax(1).XLim], [ci(1) ci(1)], '--r', [ax(1).XLim], [ci(2) ci(2)], '--r');

% Display grouped barplots of harmaline days
ax(2) = subplot(1,3,2:3);
set(ax(2), 'XLim', [0.5143 4.4857]);
plot(ax(2), [ax(2).XLim], [ci(1) ci(1)], '--r', [ax(2).XLim], [ci(2) ci(2)], '--r');
hold on
bar(har); legend(ax(2).Children(1:2), 'expos 2', 'expos 1', ...
    'Location', 'northeastoutside');
title('Harmaline exposure TS task days');
ax(2).XTickLabel = {'2 mg/kg'; ''; '4 mg/kg'; ''; '6 mg/kg'; ''; '8 mg/kg'};
grid on;

set(ax(:), 'YLim', [0 Ymax]);



%% % of All reach attempts that actually succeeded in being a reach

fieldStr = 'x_invalidInitVSallReachAttempts';

har = [T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 2) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 1) & (T.harDelivered == true)), ... 
    T.(fieldStr)((T.harDose == 4) & (T.doseExpos == 2) & (T.harDelivered == true)); ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 6) & (T.doseExpos == 2) & (T.harDelivered == true)); ... 
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 1) & (T.harDelivered == true)), ...
    T.(fieldStr)((T.harDose == 8) & (T.doseExpos == 2) & (T.harDelivered == true));];
har = 100 - har;

m = T.(fieldStr)(T.harDelivered == false); % fill this manually with just maintenance day successTrial data
m = 100 - m;

figure;
% Display violin plot of maintenance days
ax(1) = subplot(1,3,1); 
violinplot(m, [], 'Width', 0.3);
title('Naive-maintenace TS task days'); 
% ylabel('# target-reaches'); grid on;

% Display grouped barplots of harmaline days
ax(2) = subplot(1,3,2:3);
bar(har); legend('expos 1', 'expos 2');
title('Harmaline exposure TS task days');
ax(2).XTickLabel = {'2 mg/kg'; '4 mg/kg'; '6 mg/kg'; '8 mg/kg'};
grid on;

set(ax(:), 'YLim', [0 100]);
set(gcf, 'Position', [2008 360 842 420]);



