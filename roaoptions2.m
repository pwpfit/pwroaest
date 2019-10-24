classdef roaoptions2 < roaoptions
% Extended options for polynomial ROA estimation.

properties
    log = 'none';
    logpath;
end

methods 
    function opt = roaoptions2(varargin)
        opt@roaoptions(varargin{:});
        
        if isempty(opt.logpath) && ~strcmp(opt.log,'none')
            opt.logpath = 'data/';
        end
    end
    
    % Set: log
    function opt = set.log(opt,value)
        AllowableVal = {'none' 'step' 'result'};
        if ischar(value) && any(strcmp(value,AllowableVal))
            opt.log = value;
        else
            error('log must be one of ''none,'' ''step,'' or ''result''.');
        end
    end
    
    % Set: logpath
    function opt = set.logpath(opt,value)
        if ischar(value)
            opt.logpath = {value};
        elseif iscell(value) && ~isempty(value) && ischar(value{1})
            opt.logpath = value;
        else
            error('log path must be character array.');
        end
    end 
end

end
    