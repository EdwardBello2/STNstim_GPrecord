% function for up-sampling non-uniformly sampled time-series data into
% uniformly-sampled data at a higher sampling frequency. 
%
%
% Syntax:
% [sigUniform, tUniform] = upsampUniform(sigIrreg, tIrreg, fsDesired)
%
%
% Description:
% upsampUniform returns resampled data and resampled time array, according
% to the desired sampling frequency "fsDesired". Note that the input
% irregular time array should be entered in seconds. Aliasing will not be a 
% problem and edge-effects are taken care of. 
%

% Author: Ed Bello
% Last Updated: 2019/08/14
%
% NOTE: I didn't write this with down-sampling in mind, so can't guarantee
% that things will work well for downsampling/filtering, but should work in
% principle...

function [sigUniform, tUniform] = upsampUniform(sigIrreg, tIrreg, fsDesired)
% Resample it to be uniformly-sampled
% while specifying P and Q for an optimal intermediate interpolation grid.
% Since we're upsampling, aliasing is not a problem, and furthermore
% "resample.m" default anti-alising filter is run, but doesn't have much of
% an effect. Edge-effects are controlled by detrending the data,
% resampling, and then adding back in the trend. 

% % Make sure that sigIrreg and tIrreg are both column vectors
% if ~isvector(sigIrreg) | ~isvector(tIrreg)
%     error('inputs must be vectors!')
%     
% end

if size(sigIrreg, 1) < size(sigIrreg, 2), sigIrreg = sigIrreg'; end
if size(tIrreg, 1) < size(tIrreg, 2), tIrreg = tIrreg'; end

% accIrreg = file.scalars.UDP1.data(1,:);
% tIrreg = file.scalars.UDP1.ts;

% get the average sampling frequency based on the irregular sampling times:
fsNominal = (numel(tIrreg) - 1) / (tIrreg(end) - tIrreg(1));


% Tighter interpolation-grid with endpoint-effect correction:
[p, q] = rat(fsNominal / fsDesired);


nCols = size(sigIrreg, 2);
for iCol = 1:nCols
    % compute slope and offset (y = a1 x + a2)
    c(1) = (sigIrreg(end,iCol)-sigIrreg(1,iCol)) ./ (tIrreg(end)-tIrreg(1));
    c(2) = sigIrreg(1,iCol);

    % detrend the signal
    % sigDetrend = sigIrreg - polyval(c, tIrreg);
    % 
    % y = (c(1) * tIrreg) - c(2);
    sigIrregDetrend(:,iCol) = sigIrreg(:,iCol) - polyval(c, tIrreg);

    % % Resample and add back in the trend
    [sigRegPQdetr(:,iCol),tUniform] = resample(sigIrregDetrend(:,iCol), tIrreg, fsDesired, p, q);
    % 
    sigUniform(:,iCol) = sigRegPQdetr(:,iCol) + polyval(c, tUniform);
    
end

% % compute slope and offset (y = a1 x + a2)
% c(1,:) = (sigIrreg(end,:)-sigIrreg(1,:)) ./ (tIrreg(end)-tIrreg(1));
% c(2,:) = sigIrreg(1);
% 
% % detrend the signal
% % sigDetrend = sigIrreg - polyval(c, tIrreg);
% % 
% % y = (c(1) * tIrreg) - c(2);
% sigIrregDetrend = sigIrreg - polyval(c, tIrreg);
% 
% % % Resample and add back in the trend
% [sigRegPQdetr,tUniform] = resample(sigIrregDetrend, tIrreg, fsDesired, p, q);
% % 
% sigUniform = sigRegPQdetr + polyval(c, tUniform);

% 
% % Resample and add back in the trend
% [sigRegPQdetr,tUniform] = resample(sigIrreg, tIrreg, fsDesired, p, q);
% 
% sigUniform = sigRegPQdetr + polyval(c, tUniform);
% 


end