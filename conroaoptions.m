classdef conroaoptions < roaoptions
% Options for ROA estimation under constraint control.
%
%% Usage & description
%
%   opts = conroaoptions(f, g, x, ...)
%
% with inputs
%       -f: polynomial vector field
%       -x: state-space vector as PVAR
%       -u: input vector as PVAR
%
%   Name, Value:
%       -c: constraint polynomial, vector field
%       -zg: Monomials for constraint multiplier; zg is a column vector of
%            monomials used to specifiy sg(x) in the Gram matrix form.
%            [default = monomials(x, 1:2)]
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-12-01
% * Changed:    2018-12-01
%
%% See also
%
% See ROAOPTIONS
%%

properties
    %f -- inherited from ROAOPTIONS
    %x -- inherited from ROAOPTIONS
    c;
    u;
    zg;
    
    zK;
    Kin;
end

methods
    function opt = conroaoptions(f, x, u, varargin)
        if nargin < 4 || ischar(u)
            varargin = [{u} varargin];
            u = [];
        end
        
        opt@roaoptions(f, x, varargin{:});
        
        opt.u = u;

        if isempty(opt.zg)
            opt.zg = monomials(opt.x, 0:2);
        end
    end
    
    function opt = set.c(opt,value)
        if isa(c,'polynomial')
            opt.c = value;
        end
    
    function opt = set.zg(opt,value)
        if ismonom(value)
            opt.zg = value;
        else
            error('zg must be vector of monomials.');
        end
    end
end

end