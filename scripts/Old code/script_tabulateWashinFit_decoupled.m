% script to put cftool fitted results into easy tabular format

order = 8;
dose = 8;
exposure = 2;

% load the cftool session and export the variables to workspace
SSE = goodness.sse;
Rsquare = goodness.rsquare;
RsquareAdj = goodness.adjrsquare;
RMSE = goodness.rmse;

%
coeff = coeffvalues(fittedmodel);
CI = confint(fittedmodel);

a = coeff(1);
alower = CI(1,1);
aupper = CI(2,1);

b = coeff(2);
blower = CI(1,2);
bupper = CI(2,2);

k = coeff(3);
klower = CI(1,3);
kupper = CI(2,3);



T = [table(order), table(dose), table(exposure), ...
    table(a), table(alower), table(aupper), ...
    table(b), table(blower), table(bupper), ...
    table(k), table(klower), table(kupper), ...
    table(SSE), table(Rsquare), table(RsquareAdj), table(RMSE)]

open T
clear fittedmodel goodness output