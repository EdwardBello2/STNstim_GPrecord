function [sigFilt] = filtLineNoise(sig,fs,nHarmonics)
% notch-filter 60 hz line noise, and any harmonics



fundamental = 60; % hz

% if user does not specify nHarmonics, notch only fundamental freq
if nargin < 3 
    nHarmonics = 0;
    
end

assert((fundamental * (nHarmonics+1)) <= (fs / 2), ...
    'Too many harmonics, will not fit within Nyquist frequency');

sigProc = sig;
bw = 4 / (fs / 2);
for i = 1:(nHarmonics + 1)
    wo = (i * fundamental) / (fs / 2);  
    [b,a] = iirnotch(wo, bw);
    sigProc = filtfilt(b, a, sigProc);
    
end

sigFilt = sigProc;

end