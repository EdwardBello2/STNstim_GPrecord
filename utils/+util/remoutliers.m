function [B, TF] = remoutliers(A, varargin)
% Uses Median Absolute Deviation for outliers, anything more than 3*MAD is
% considered outlier. I made this util function in pathetic mockery of matlab's
% official rmoutliers, but alas I do not currently have 2018b...
% Ed Bello


p = inputParser;

addRequired(p, 'A');
addParameter(p, 'bound', 'both', @ischar);
addParameter(p, 'MADthresh', 3, @isnumeric);

parse(p, A, varargin{:});

% A = p.Results.A;
bound = p.Results.bound;
MADthresh = p.Results.MADthresh;

if size(A, 1) < size(A, 2) % assert column vector
    A = A';
    
end

c = -1/(sqrt(2)*erfcinv(3/2));

MAD = c * median(abs(A - median(A, 'omitnan')), 'omitnan');

switch bound
    case 'both'
        TF = (A < (median(A, 'omitnan') - MADthresh*MAD)) | (A > (median(A, 'omitnan') + MADthresh*MAD));
        
    case 'upper'
        TF = (A > (median(A, 'omitnan') + MADthresh*MAD));
        
    case 'lower'
        TF = (A < (median(A, 'omitnan') - MADthresh*MAD));
        
    otherwise
        TF = (A < (median(A, 'omitnan') - MADthresh*MAD)) | (A > (median(A, 'omitnan') + MADthresh*MAD));
        
end

B = A;
B(TF) = [];

end