function outVect = ArtFilt(inVect,samp,HPass,LPass, isNotch)

nPoles = 4;

if (HPass && LPass)
    Wn=[HPass LPass]/(samp/2);
    [b,a]=butter(nPoles, Wn, 'bandpass');
elseif (HPass)
    Wn=[HPass]/(samp/2);
    [b,a]=butter(nPoles, Wn,'high');
elseif (LPass)
    Wn=[LPass]/(samp/2);
    [b,a]=butter(nPoles, Wn,'low');
end

if (HPass || LPass)
    inVect=filtfilt(b,a,inVect);
end

if (isNotch)
    wo = 50/(samp/2);  bw = wo/35;
    [bnotch,anotch] = iirnotch(wo,bw);
    outVect=filtfilt(bnotch,anotch,inVect);
else
    outVect=inVect;
end