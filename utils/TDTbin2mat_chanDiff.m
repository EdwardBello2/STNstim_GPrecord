% Load tdt data identical to the usual TDTbin2mat in all respects except that it takes
% the differential signal between two channels specified. 
%
%
% SYNTAX:
% tdt = TDTbin2mat_chanDiff(blockpath, chPair)
% tdt = TDTbin2mat_chanDiff(__, Name, Value)
%
% (see TDTbin2mat help for Name-Value pairs)
%
% 
% DESCRIPTION:
% Simply specify the full path to your TDT block in blockpath, and specify
% the two channels for subtracting in chPair (smaller channel index first),
% followed by Name-Value pair arguments according to standard calls of
% TDTbin2mat. 
%
%
% INPUTS:
%   blockpath - fullpath to directory holding all TDT blocks, char array
%      chPair - n x 1 array of integer indices representing channels to
%               subtract
%
%
% OUTPUTS: 
% tdt - struct with one channel of stream data, so that data can be
%       accessed with dot-notation (e.g. tdtConcat.streams.data).
%
% 
%
% EXAMPLES:
% 


% Author: Ed Bello
% Created: 3/26/2020
%
% 
% ChangeLog:

% 
%
% TO-DO


function  tdt2 = TDTbin2mat_chanDiff(blockpath,chPair, varargin)



%% Initial section identical to TDTbin2mat code (thanks TDT!)

% defaults
BITWISE  = '';
CHANNEL  = 0;
COMBINE  = {};
HEADERS  = 0;
NODATA   = false;
RANGES   = [];
STORE    = '';
T1       = 0;
T2       = 0;
TYPE     = 1:5;
VERBOSE  = 0;
SORTNAME = 'TankSort';

% VALID_PARS = {'BITWISE','CHANNEL','HEADERS','NODATA','RANGES','STORE', ...
%     'T1','T2','TYPE','VERBOSE','SORTNAME','COMBINE'};
VALID_PARS = {'BITWISE','HEADERS','NODATA','RANGES','STORE', ...
    'T1','T2','TYPE','VERBOSE','SORTNAME','COMBINE'};

% parse varargin
for ii = 1:2:length(varargin)
    if ~ismember(upper(varargin{ii}), VALID_PARS)
        error('%s is not a valid parameter. See help TDTbin2mat.', upper(varargin{ii}));
    end
    eval([upper(varargin{ii}) '=varargin{ii+1};']);
end



%% NRTL code

% make sure that chPair is a two-element array
assert((numel(chPair) == 2) & (isvector(chPair)), ...
    'chPair must be a two-element vector-array'); 

tdt1 = TDTbin2mat(blockpath, 'CHANNEL', chPair(1), ...
    'COMBINE', COMBINE, ...
    'HEADERS', HEADERS, ...
    'NODATA', NODATA, ...
    'RANGES', RANGES, ...
    'STORE', STORE, ...
    'T1', T1, ...
    'T2', T2, ...
    'TYPE', TYPE, ...
    'VERBOSE', VERBOSE, ...
    'SORTNAME', SORTNAME);

tdt2 = TDTbin2mat(blockpath, 'CHANNEL', chPair(2), ...
    'COMBINE', COMBINE, ...
    'HEADERS', HEADERS, ...
    'NODATA', NODATA, ...
    'RANGES', RANGES, ...
    'STORE', STORE, ...
    'T1', T1, ...
    'T2', T2, ...
    'TYPE', TYPE, ...
    'VERBOSE', VERBOSE, ...
    'SORTNAME', SORTNAME);

tdt2.streams.(STORE).channel = [tdt1.streams.(STORE).channel, ...
    tdt2.streams.(STORE).channel];

tdt2.streams.(STORE).data = tdt2.streams.(STORE).data - ...
    tdt1.streams.(STORE).data;




end