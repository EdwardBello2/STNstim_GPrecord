function [stim_est,filt_sig,runTime] = aFIR_filtfilt(raw_data,stimTime,order,LR)
%% Adaptive FIR filter for electrical stimulation artifact removal
% Authors: Kenneth Louie[1], Vivek Nagaraj[2], & Edward Bello[1]
% Affliation:   [1] University of Minnesota Department of Biomedical
%                   Engineering
%               [2] University of Minnesota Department of Neuroscience
% Version: 1.0  Initial
%          1.1  Added run time output
%          1.2  Added in filtfilt-like functionality
% Date Completed: v1.0 10/25/2016
%                 v1.1 12/12/2016
%                 v1.2 05/19/2020
% DESCRIPTION:
%  Uses the least mean squares (LMS) algorithm to tune FIR filter
%  coefficients. This is done to minimize the error between the desired and
%  actual signal.
% INPUTS:   raw_data    [=] Output signal (un processed data), row vector
%                           or column vector
%           stimTime    [=] Same length as raw_data. This is a vector of
%                           0's and 1's, indicating a stimulation was
%                           delivered
%           order       [=] Order of the FIR filter
%           LR          [=] The learning rate for the FIR filter
%                           coefficients
% OUTPUTS:  stim_est    [=] Estimated stimulus model. This is what is
%                           subtracted from the raw data.
%           filt_sig    [=] The filtered signal after the FIR filter was
%                           applied
%           runTime     [=] Time it took to run the filter on the data

% Obtain reverse version of stimTimes, with 1's on other end of filter 
stimTimeRev = zeros(1, length(stimTime));
stimIdx = find(stimTime);
stimTimeRev(stimIdx+order-1) = 1;

% First flip the order of input data, run this filter backwards on first
% pass
raw_dataRev = flip(raw_data);
stimTimeRev = flip(stimTimeRev);


% Miscellaneous variables
N = length(raw_data);   % Length of the raw data
b = zeros(order, 1);     % Initial FIR filter coefficients

% Output vectors
stim_est = zeros(N, 1);
filt_sig = zeros(N, 1);

% Run filter backwards first
A = tic;
[~, ~, b] = aFIRfilt(raw_dataRev, stimTimeRev, order, LR, ...
    b, N, stim_est, filt_sig);


% Now flip the FIR filter coefficients and run it forwards, reset outputs
% to zeros
b = flip(b);
stim_est = zeros(N, 1);
filt_sig = zeros(N, 1);

[stim_est, filt_sig, ~] = aFIRfilt(raw_data, stimTime, order, LR, ...
    b, N, stim_est, filt_sig);

runTime = toc(A);
end

function [stim_est,filt_sig,b] = aFIRfilt(raw_data,stimTime,order,LR, b, N, stim_est, filt_sig)

for i=order+1:N
    % B Coefficients
    for j=1:order
        stim_est(i) = stim_est(i) + stimTime(i-j+1)*b(j);
    end
    
    % Subtract the estimated stimulus model from the raw data
    filt_sig(i) = raw_data(i) - stim_est(i);
    
    % Update coefficients
    for j=1:order
        b(j) = b(j) + LR*filt_sig(i)*stimTime(i-j+1);
    end
end

end