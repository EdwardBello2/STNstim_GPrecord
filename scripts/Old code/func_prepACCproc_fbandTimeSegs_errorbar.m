function [tmid,pow,cineg,cipos,tneg,tpos] = func_prepACCproc_fbandTimeSegs_errorbar(tremPow, tstSeg, cfg, statFunLabel)
%% display chunks of recorded tremor power in time time showing summary 
% statistic of tremor power and confidence interval of that statistic for
% a given chunk of data, using errorbar function

if exist('statFunLabel', 'var')
    switch statFunLabel
        case 'mean'
            myStatFun = @(x)mean(x); 
            
        case 'median'
            myStatFun = @(x)median(x); 
            
        otherwise
            myStatFun = @(x)median(x);
            
    end
  
else
    myStatFun = @(x)median(x); 
    
end

% CItype: method in bootci to use for bootstrapped estimate of confidence
% intervals
CItype = 'bca';

%%

% Get time-span of segments organized
nSegs = length(tstSeg);
for iSeg = 1:nSegs
    % get time span for each segment, time midpoint
    iTst = tstSeg{iSeg};
    tspan(iSeg) = iTst(end) - iTst(1);
    tmid(iSeg) = (iTst(end) + iTst(1)) / 2;
    
end

% convert from seconds to minutes
tmid = tmid / 60;
tspan = tspan / 60;


% Get tremor power summary statistics
for iSeg = 1:nSegs
    % distribution of tremor power for iSeg
    
    % get mean of each segment
    tremPowav(iSeg) = myStatFun(tremPow{iSeg,1});
    
    try
        % get confidence interval of that mean
        ci(1:2,iSeg) = bootci(cfg.NBOOTS, {myStatFun, tremPow{iSeg,1}}, 'type', CItype);
    
    catch ME
        disp(['Loop iteration: ' num2str(iSeg)])
        error(ME.message)
        
    end
        
    
end


    
for iSeg = 1:nSegs
%     x(iSeg) = tmid(iSeg);
    tneg(iSeg) = 0.5*tspan(iSeg);
    tpos(iSeg) = 0.5*tspan(iSeg);
    
end

pow = tremPowav;
cineg = pow - ci(1,:);
cipos = ci(2,:) - pow;


end