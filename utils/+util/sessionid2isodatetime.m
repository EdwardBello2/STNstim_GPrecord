function [isoStr] = sessionid2isodatetime(sessidStr)
% Extracts date and time info from my standard sessionID filenaming format
% (based on how TDT blocks are named) and outputs an ISO 8601
% extended formatted string for all that info. 

isoStr = ['20', sessidStr(end-12:end-11), '-', sessidStr(end-10:end-9), '-', sessidStr(end-8:end-7), ...
    'T', sessidStr(end-5:end-4), ':', sessidStr(end-3:end-2), ':', sessidStr(end-1:end), '.000Z'];

end