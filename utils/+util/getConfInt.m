function CI = getConfInt(x, CIperc)
% A little function based on Star Stirders' solution for getting confidence
% intervals on a given dataset https://www.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval



% x is a vector, matrix, or any numeric array of data. NaNs are ignored.
% p is the confidence level (ie, 95 for 95% CI)
% The output is 1x2 vector showing the [lower,upper] interval values.

% x = pearsrnd(0,1,1,4,100,1); 
% figure; histogram(x);  

p = CIperc;
% 
% CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),sum(~isnan(x(:)))-1) + mean(x(:),'omitnan'); 
% CI = CIFcn(x,95); 
% arrayfun(@(x)xline(x,'-k','tinv'),CI);

CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
CI = CIFcn(x,p); 




end