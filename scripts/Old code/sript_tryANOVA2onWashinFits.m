% script to test out simple 2-way ANOVA on harmaline washin fit params

reps = 1; % for both factors (exposure, dose) I only have 1 example for each combo ([1,2], [2,4,6,8])

% Fill y manually with the fit params
% Rows are dose (2,4,6,8)
% Cols are Exposure (1,2)

% High-beta "k" coefficient (final asymptote)
y = [1229, 378.3; 260.2, 378.3; 416.4, 630.5; 433, 220.5]
[p, tbl, stats] = anova2(y, reps)



% High-beta "v" coefficient (final asymptote)
y = [3.898, 4.458; 5.615, 6.184; 4.256, 5.483; 4.227, 5.937]
[p, tbl, stats] = anova2(y, reps)



% Low-beta "k" coefficient (final asymptote)
y = [10, 10; 180.5, 571.6; 1546, 920.4; 389.3, 971.9]
[p, tbl, stats] = anova2(y, reps)



% Low-beta "v" coefficient (final asymptote)
y = [1.67, 2.896; 3.516, 4.294; 2.153, 4.993; 3.616, 4.813]
[p, tbl, stats] = anova2(y, reps)