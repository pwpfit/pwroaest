classdef droaoptions < roaoptions2
% Options for discrete ROA estimation.
%
%% Usage & description
%
%   opts = droaoptions(f, x, ...)
%
% where f is the continuous system's equation s.t.
%
%   x+ = x + tau*f(x)
%
% with inputs
%       -f, x: see ROAOPTIONS
%
%   Name, Value:
%       -tau: sample time [default = 0]
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-05-22
% * Changed:    2018-10-28
%
%% See also
%
% See ROAOPTIONS
%%

properties
    %f -- inherited from ROAOPTIONS
    %x -- inherited form ROAOPTIONS
    tau;
end

methods
    function opt = droaoptions(varargin)
        opt@roaoptions2(varargin{:});
        
        if isempty(opt.tau)
            opt.tau = 0;
        end
    end
    
    function opt = set.tau(opt,value)
        if isnumeric(value)
            opt.tau = value;
        else
            error('Sample time must be a number.');
        end
    end
end

end