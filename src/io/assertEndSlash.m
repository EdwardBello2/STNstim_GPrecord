function str = assertEndSlash(str)
if ~strcmp(str(end), '\'), str = [str '\']; end

end